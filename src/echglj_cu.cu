#include "add.h"
#include "box.h"
#include "couple.h"
#include "echglj.h"
#include "glob.spatial.h"
#include "image.h"
#include "launch.h"
#include "mathfunc.h"
#include "md.h"
#include "osrw.h"
#include "pmestuf.h"
#include "seq_switch.h"
#include "switch.h"
#include "tool/gpu_card.h"


namespace tinker {
namespace {
cudaEvent_t echglj_start;
cudaStream_t echglj_stream;
}


void echglj_cu_data(rc_op op)
{
   if (op bitand rc_dealloc) {
      cudaEventDestroy(echglj_start);
      cudaEventDestroy(echglj_event);
      echglj_event = nullptr;
      use_echglj_event = false;
      echglj_stream = nullptr;
   }


   if (op bitand rc_alloc) {
      cudaEventCreateWithFlags(&echglj_start, cudaEventDisableTiming);
      cudaEventCreateWithFlags(&echglj_event, cudaEventDisableTiming);
      use_echglj_event = true;
      echglj_stream = nullptr; // the default 0 CUDA stream
   }
}
}


extern "C"
{
   using real = tinker::real;


   struct RAD_GEOM
   {
      __device__
      static real avg(real ra, real rb)
      {
         return 2 * REAL_SQRT(ra * rb);
      }
   };


   struct RAD_ARITH
   {
      __device__
      static real avg(real ra, real rb)
      {
         return (ra + rb);
      }
   };


   struct EPS_GEOM
   {
      __device__
      static real avg(real a, real b)
      {
         return REAL_SQRT(a * b);
      }
   };


   struct EPS_ARITH
   {
      __device__
      static real avg(real a, real b)
      {
         return 0.5f * (a + b);
      }
   };
}


namespace tinker {
template <bool DO_G, class ETYP, class RADRULE, class EPSRULE>
__device__
void pair_chglj(real r, real r2, real& restrict invr, //
                real cscale, real chgi, real chgk, real ebuffer, real f,
                real aewald, real eccut, real ecoff, real& restrict ec,
                real& restrict dec, //
                real vscale, real radi, real epsi, real radk, real epsk,
                real evcut, real evoff, real& restrict ev, real& restrict dev)
{
   if CONSTEXPR (DO_G)
      invr = REAL_RECIP(r);


   // vdw
   ev = 0;
   dev = 0;
   //*
   if (r <= evoff) {
      real eps = EPSRULE::avg(epsi, epsk) * vscale;
      real rv = RADRULE::avg(radi, radk);
      real rv2_rik2 = rv * rv * REAL_RECIP(r2);
      real p6 = rv2_rik2 * rv2_rik2 * rv2_rik2;
      ev = eps * p6 * (p6 - 2);
      if CONSTEXPR (DO_G) {
         dev = eps * p6 * (p6 - 1) * (-12 * invr);
      }
      // taper
      if (r > evcut) {
         real taper, dtaper;
         switch_taper5<DO_G>(r, evcut, evoff, taper, dtaper);
         if CONSTEXPR (DO_G)
            dev = ev * dtaper + dev * taper;
         ev = ev * taper;
      }
   }
   // */


   // charge
   ec = 0;
   dec = 0;
   //*
   constexpr bool taper_flag = eq<ETYP, NON_EWALD_TAPER>();
   if (r <= ecoff) {
      if CONSTEXPR (eq<ETYP, EWALD>()) {
         real fik = f * chgi * chgk;
         real rew = aewald * r;
         real erfterm = REAL_ERFC(rew);
         real rb = r + ebuffer;
         real invrb = REAL_RECIP(rb);
         ec = fik * invrb * erfterm;
         if CONSTEXPR (DO_G) {
            real invrb2 = invrb * invrb;
            dec = erfterm * invrb2;
            dec += (2 * aewald / sqrtpi) * REAL_EXP(-rew * rew) * invr;
            dec *= -fik;
         }
      } else if CONSTEXPR (eq<ETYP, NON_EWALD>() || taper_flag) {
         real fik = cscale * f * chgi * chgk;
         real rb = r + ebuffer;
         real invrb = REAL_RECIP(rb);
         ec = fik * invrb;
         if CONSTEXPR (DO_G) {
            real invrb2 = invrb * invrb;
            dec = -fik * invrb2;
         }


         if CONSTEXPR (taper_flag) {
            // shift energy
            real shift = fik * 2 * REAL_RECIP(eccut + ecoff);
            ec -= shift;
            if (r > eccut) {
               real taper, dtaper;
               switch_taper5<DO_G>(r, eccut, ecoff, taper, dtaper);


               real trans, dtrans;
               real coef = fik * (REAL_RECIP(eccut) - REAL_RECIP(ecoff)) *
                  REAL_RECIP((real)9.3);
               real invs = REAL_RECIP(ecoff - eccut);
               real x = (r - eccut) * invs;
               real y = (1 - x) * x;
               trans = coef * y * y * y * (64 - 25 * x);
               if CONSTEXPR (DO_G) {
                  dtrans = y * y * (25 * x - 12) * (7 * x - 16);
                  dtrans *= coef * invs;
                  dec = ec * dtaper + dec * taper + dtrans;
               }
               ec = ec * taper + trans;
            }
         }
      }
   } // end if (include ec)
   // */
}


//====================================================================//


#define ECHGLJPARAS                                                            \
   size_t bufsize, energy_buffer restrict ec, virial_buffer restrict vir_ec,   \
      grad_prec *restrict decx, grad_prec *restrict decy,                      \
      grad_prec *restrict decz, real eccut, real ecoff, real ebuffer, real f,  \
      const real *restrict pchg, energy_buffer restrict ev,                    \
      virial_buffer restrict vir_ev, grad_prec *restrict devx,                 \
      grad_prec *restrict devy, grad_prec *restrict devz, real *restrict rad,  \
      real *restrict eps, real evcut, real evoff, TINKER_IMAGE_PARAMS


template <class Ver, class ETYP, class RADRULE, class EPSRULE>
__global__
void echglj_cu1(ECHGLJPARAS, const Spatial::SortedAtom* restrict sorted,
                int niak, const int* restrict iak, const int* restrict lst,
                int n, real aewald, const int (*restrict i12)[couple_maxn12])
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);
   const int offset = ithread & (bufsize - 1);


   struct Data
   {
      real x, y, z, c, rad, eps;
      grad_prec fcx, fcy, fcz, fvx, fvy, fvz;
   };
   // real == float, grad_prec == float, NByte: 4 * 12 = 48
   // real == float, grad_prec == fixed, NByte: 4 * 6 + 8 * 6 = 72
   // real == double, grad_prec == double, NByte: 8 * 12 = 96
   // real == double, grad_prec == fixed, NByte: 8 * 12 = 96
   __shared__ Data data[BLOCK_DIM];


   // thread local variables
   using ebuf_prec = energy_buffer_traits::type;
   using vbuf_prec = virial_buffer_traits::type;
   MAYBE_UNUSED ebuf_prec ectl, evtl;
   MAYBE_UNUSED vbuf_prec vctlxx, vctlyx, vctlzx, vctlyy, vctlzy, vctlzz;
   MAYBE_UNUSED vbuf_prec vvtlxx, vvtlyx, vvtlzx, vvtlyy, vvtlzy, vvtlzz;
   Data idat;


   for (int iw = iwarp; iw < niak; iw += nwarp) {
      if CONSTEXPR (do_e) {
         ectl = 0;
         evtl = 0;
      }
      if CONSTEXPR (do_g) {
         idat.fcx = 0;
         idat.fcy = 0;
         idat.fcz = 0;
         idat.fvx = 0;
         idat.fvy = 0;
         idat.fvz = 0;
         data[threadIdx.x].fcx = 0;
         data[threadIdx.x].fcy = 0;
         data[threadIdx.x].fcz = 0;
         data[threadIdx.x].fvx = 0;
         data[threadIdx.x].fvy = 0;
         data[threadIdx.x].fvz = 0;
      }
      if CONSTEXPR (do_v) {
         vctlxx = 0;
         vctlyx = 0;
         vctlzx = 0;
         vctlyy = 0;
         vctlzy = 0;
         vctlzz = 0;
         vvtlxx = 0;
         vvtlyx = 0;
         vvtlzx = 0;
         vvtlyy = 0;
         vvtlzy = 0;
         vvtlzz = 0;
      }


      int atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
      idat.x = sorted[atomi].x;
      idat.y = sorted[atomi].y;
      idat.z = sorted[atomi].z;
      int i = sorted[atomi].unsorted;
      idat.c = pchg[i];
      idat.rad = rad[i];
      idat.eps = eps[i];


      int shatomk = lst[iw * WARP_SIZE + ilane];
      data[threadIdx.x].x = sorted[shatomk].x;
      data[threadIdx.x].y = sorted[shatomk].y;
      data[threadIdx.x].z = sorted[shatomk].z;
      int shk = sorted[shatomk].unsorted;
      data[threadIdx.x].c = pchg[shk];
      data[threadIdx.x].rad = rad[shk];
      data[threadIdx.x].eps = eps[shk];


      int cpli[couple_maxn12];
      if (i12) {
         #pragma unroll
         for (int ic = 0; ic < couple_maxn12; ++ic) {
            cpli[ic] = i12[i][ic];
         }
      }


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
         int k = __shfl_sync(ALL_LANES, shk, srclane);
         real xr = idat.x - data[klane].x;
         real yr = idat.y - data[klane].y;
         real zr = idat.z - data[klane].z;


         real vscale = 1;
         bool ik_bond = false;
         if (i12) {
            #pragma unroll
            for (int ic = 0; ic < couple_maxn12; ++ic) {
               ik_bond = ik_bond || (cpli[ic] == k);
            }
         }
         if (ik_bond)
            vscale = 0;


         MAYBE_UNUSED real dedxc = 0, dedyc = 0, dedzc = 0;
         MAYBE_UNUSED real dedxv = 0, dedyv = 0, dedzv = 0;


         real r2 = image2(xr, yr, zr);
         if (atomi < atomk) {
            MAYBE_UNUSED real invr;
            real r = REAL_SQRT(r2);
            real ecik, decik;
            real evik, devik;
            pair_chglj<do_g, ETYP, RADRULE, EPSRULE>(
               r, r2, invr, //
               1, idat.c, data[klane].c, ebuffer, f, aewald, eccut, ecoff, ecik,
               decik, //
               vscale, idat.rad, idat.eps, data[klane].rad, data[klane].eps,
               evcut, evoff, evik, devik);
            if CONSTEXPR (do_e) {
               ectl += cvt_to<ebuf_prec>(ecik);
               evtl += cvt_to<ebuf_prec>(evik);
            }
            if CONSTEXPR (do_g) {
               decik *= invr;
               dedxc = decik * xr;
               dedyc = decik * yr;
               dedzc = decik * zr;
               devik *= invr;
               dedxv = devik * xr;
               dedyv = devik * yr;
               dedzv = devik * zr;
               if CONSTEXPR (do_v) {
                  vctlxx += cvt_to<vbuf_prec>(xr * dedxc);
                  vctlyx += cvt_to<vbuf_prec>(yr * dedxc);
                  vctlzx += cvt_to<vbuf_prec>(zr * dedxc);
                  vctlyy += cvt_to<vbuf_prec>(yr * dedyc);
                  vctlzy += cvt_to<vbuf_prec>(zr * dedyc);
                  vctlzz += cvt_to<vbuf_prec>(zr * dedzc);
                  vvtlxx += cvt_to<vbuf_prec>(xr * dedxv);
                  vvtlyx += cvt_to<vbuf_prec>(yr * dedxv);
                  vvtlzx += cvt_to<vbuf_prec>(zr * dedxv);
                  vvtlyy += cvt_to<vbuf_prec>(yr * dedyv);
                  vvtlzy += cvt_to<vbuf_prec>(zr * dedyv);
                  vvtlzz += cvt_to<vbuf_prec>(zr * dedzv);
               }
            }
         } // if (include)


         if CONSTEXPR (do_g) {
            idat.fcx += cvt_to<grad_prec>(dedxc);
            idat.fcy += cvt_to<grad_prec>(dedyc);
            idat.fcz += cvt_to<grad_prec>(dedzc);
            data[klane].fcx -= cvt_to<grad_prec>(dedxc);
            data[klane].fcy -= cvt_to<grad_prec>(dedyc);
            data[klane].fcz -= cvt_to<grad_prec>(dedzc);
            idat.fvx += cvt_to<grad_prec>(dedxv);
            idat.fvy += cvt_to<grad_prec>(dedyv);
            idat.fvz += cvt_to<grad_prec>(dedzv);
            data[klane].fvx -= cvt_to<grad_prec>(dedxv);
            data[klane].fvy -= cvt_to<grad_prec>(dedyv);
            data[klane].fvz -= cvt_to<grad_prec>(dedzv);
         }
      }


      if CONSTEXPR (do_e) {
         atomic_add(ectl, ec, offset);
         atomic_add(evtl, ev, offset);
      }
      if CONSTEXPR (do_g) {
         atomic_add(idat.fcx, decx, i);
         atomic_add(idat.fcy, decy, i);
         atomic_add(idat.fcz, decz, i);
         atomic_add(data[threadIdx.x].fcx, decx, shk);
         atomic_add(data[threadIdx.x].fcy, decy, shk);
         atomic_add(data[threadIdx.x].fcz, decz, shk);
         atomic_add(idat.fvx, devx, i);
         atomic_add(idat.fvy, devy, i);
         atomic_add(idat.fvz, devz, i);
         atomic_add(data[threadIdx.x].fvx, devx, shk);
         atomic_add(data[threadIdx.x].fvy, devy, shk);
         atomic_add(data[threadIdx.x].fvz, devz, shk);
      }
      if CONSTEXPR (do_v) {
         atomic_add(vctlxx, vctlyx, vctlzx, vctlyy, vctlzy, vctlzz, vir_ec,
                    offset);
         atomic_add(vvtlxx, vvtlyx, vvtlzx, vvtlyy, vvtlzy, vvtlzz, vir_ev,
                    offset);
      }
   } // end for (iw)
}


template <class Ver, class ETYP, class RADRULE, class EPSRULE>
__global__
void echglj_cu2(ECHGLJPARAS, const real* restrict x, const real* restrict y,
                const real* restrict z, int ncvexclude,
                const int (*restrict cvexclude)[2],
                const real (*restrict cvexclude_scale)[2])
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < ncvexclude;
        ii += blockDim.x * gridDim.x) {
      int offset = ii & (bufsize - 1);


      int i = cvexclude[ii][0];
      int k = cvexclude[ii][1];
      real cscale = cvexclude_scale[ii][0];
      real vscale = cvexclude_scale[ii][1];


      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real ci = pchg[i];
      real radi = rad[i];
      real epsi = eps[i];


      real xr = xi - x[k];
      real yr = yi - y[k];
      real zr = zi - z[k];
      real ck = pchg[k];
      real radk = rad[k];
      real epsk = eps[k];


      real r2 = image2(xr, yr, zr);
      real r = REAL_SQRT(r2);


      MAYBE_UNUSED real invr, ecik, decik, evik, devik;
      pair_chglj<do_g, ETYP, RADRULE, EPSRULE>(
         r, r2, invr,                                              //
         cscale, ci, ck, ebuffer, f, 0, eccut, ecoff, ecik, decik, //
         vscale, radi, epsi, radk, epsk, evcut, evoff, evik, devik);


      if CONSTEXPR (do_e) {
         atomic_add(ecik, ec, offset);
         atomic_add(evik, ev, offset);
      }
      if CONSTEXPR (do_g) {
         decik *= invr;
         real dedxc = decik * xr;
         real dedyc = decik * yr;
         real dedzc = decik * zr;
         devik *= invr;
         real dedxv = devik * xr;
         real dedyv = devik * yr;
         real dedzv = devik * zr;
         atomic_add(dedxc, decx, i);
         atomic_add(dedyc, decy, i);
         atomic_add(dedzc, decz, i);
         atomic_add(-dedxc, decx, k);
         atomic_add(-dedyc, decy, k);
         atomic_add(-dedzc, decz, k);
         atomic_add(dedxv, devx, i);
         atomic_add(dedyv, devy, i);
         atomic_add(dedzv, devz, i);
         atomic_add(-dedxv, devx, k);
         atomic_add(-dedyv, devy, k);
         atomic_add(-dedzv, devz, k);
         if CONSTEXPR (do_v) {
            real vxxc = xr * dedxc;
            real vyxc = yr * dedxc;
            real vzxc = zr * dedxc;
            real vyyc = yr * dedyc;
            real vzyc = zr * dedyc;
            real vzzc = zr * dedzc;
            real vxxv = xr * dedxv;
            real vyxv = yr * dedxv;
            real vzxv = zr * dedxv;
            real vyyv = yr * dedyv;
            real vzyv = zr * dedyv;
            real vzzv = zr * dedzv;
            atomic_add(vxxc, vyxc, vzxc, vyyc, vzyc, vzzc, vir_ec, offset);
            atomic_add(vxxv, vyxv, vzxv, vyyv, vzyv, vzzv, vir_ev, offset);
         }
      }
   }
}


//====================================================================//


#define ECHGLJPARAS2                                                           \
   size_t bufsize, energy_buffer restrict ebuf, virial_buffer restrict vbuf,   \
      grad_prec *restrict gx, grad_prec *restrict gy, grad_prec *restrict gz,  \
      real eccut, real ecoff, real ebuffer, real f, const real *restrict pchg, \
      real *restrict rad, real *restrict eps, real evcut, real evoff,          \
      TINKER_IMAGE_PARAMS


static constexpr int CHGLJ4_BDIM = 2 * BLOCK_DIM;


template <class Ver, class ETYP, class RADRULE, class EPSRULE>
__launch_bounds__(CHGLJ4_BDIM) __global__
void echglj_cu4(ECHGLJPARAS2, const Spatial::SortedAtom* restrict sorted,
                int niak, const int* restrict iak, const int* restrict lst,
                int n, real aewald, const int (*restrict i12)[couple_maxn12])
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);
   const int offset = ithread & (bufsize - 1);


   struct Data
   {
      real x, y, z, c;
      real rad, eps;
      real fcx, fcy, fcz;
      real padding_;
   };
   __shared__ Data data[CHGLJ4_BDIM];


   // thread local variables
   using ebuf_prec = energy_buffer_traits::type;
   using vbuf_prec = virial_buffer_traits::type;
   MAYBE_UNUSED ebuf_prec ectl;
   MAYBE_UNUSED vbuf_prec vctlxx, vctlyx, vctlzx, vctlyy, vctlzy, vctlzz;
   Data idat;


   for (int iw = iwarp; iw < niak; iw += nwarp) {
      if CONSTEXPR (do_e) {
         ectl = 0;
      }
      if CONSTEXPR (do_g) {
         idat.fcx = 0;
         idat.fcy = 0;
         idat.fcz = 0;
         data[threadIdx.x].fcx = 0;
         data[threadIdx.x].fcy = 0;
         data[threadIdx.x].fcz = 0;
      }
      if CONSTEXPR (do_v) {
         vctlxx = 0;
         vctlyx = 0;
         vctlzx = 0;
         vctlyy = 0;
         vctlzy = 0;
         vctlzz = 0;
      }


      int atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
      idat.x = sorted[atomi].x;
      idat.y = sorted[atomi].y;
      idat.z = sorted[atomi].z;
      int i = sorted[atomi].unsorted;
      idat.c = pchg[i];
      idat.rad = rad[i];
      idat.eps = eps[i];


      int shatomk = lst[iw * WARP_SIZE + ilane];
      data[threadIdx.x].x = sorted[shatomk].x;
      data[threadIdx.x].y = sorted[shatomk].y;
      data[threadIdx.x].z = sorted[shatomk].z;
      int shk = sorted[shatomk].unsorted;
      data[threadIdx.x].c = pchg[shk];
      data[threadIdx.x].rad = rad[shk];
      data[threadIdx.x].eps = eps[shk];


      int cpli[couple_maxn12];
      if (i12) {
         #pragma unroll
         for (int ic = 0; ic < couple_maxn12; ++ic) {
            cpli[ic] = i12[i][ic];
         }
      }


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
         int k = __shfl_sync(ALL_LANES, shk, srclane);
         real xr = idat.x - data[klane].x;
         real yr = idat.y - data[klane].y;
         real zr = idat.z - data[klane].z;


         real vscale = 1;
         bool ik_bond = false;
         if (i12) {
            #pragma unroll
            for (int ic = 0; ic < couple_maxn12; ++ic) {
               ik_bond = ik_bond || (cpli[ic] == k);
            }
         }
         if (ik_bond)
            vscale = 0;


         MAYBE_UNUSED real dedxc = 0, dedyc = 0, dedzc = 0;


         real r2 = image2(xr, yr, zr);
         if (atomi < atomk) {
            MAYBE_UNUSED real invr;
            real r = REAL_SQRT(r2);
            real ecik, decik;
            real evik, devik;
            pair_chglj<do_g, ETYP, RADRULE, EPSRULE>(
               r, r2, invr, //
               1, idat.c, data[klane].c, ebuffer, f, aewald, eccut, ecoff, ecik,
               decik, //
               vscale, idat.rad, idat.eps, data[klane].rad, data[klane].eps,
               evcut, evoff, evik, devik);
            if CONSTEXPR (do_e) {
               ectl += cvt_to<ebuf_prec>(ecik) + cvt_to<ebuf_prec>(evik);
            }
            if CONSTEXPR (do_g) {
               decik = (decik + devik) * invr;
               dedxc = decik * xr;
               dedyc = decik * yr;
               dedzc = decik * zr;
               if CONSTEXPR (do_v) {
                  vctlxx += cvt_to<vbuf_prec>(xr * dedxc);
                  vctlyx += cvt_to<vbuf_prec>(yr * dedxc);
                  vctlzx += cvt_to<vbuf_prec>(zr * dedxc);
                  vctlyy += cvt_to<vbuf_prec>(yr * dedyc);
                  vctlzy += cvt_to<vbuf_prec>(zr * dedyc);
                  vctlzz += cvt_to<vbuf_prec>(zr * dedzc);
               }
            }
         } // end if (include)


         if CONSTEXPR (do_g) {
            idat.fcx += dedxc;
            idat.fcy += dedyc;
            idat.fcz += dedzc;
            data[klane].fcx -= dedxc;
            data[klane].fcy -= dedyc;
            data[klane].fcz -= dedzc;
         }
      }


      if CONSTEXPR (do_e) {
         atomic_add(ectl, ebuf, offset);
      }
      if CONSTEXPR (do_g) {
         atomic_add(idat.fcx, gx, i);
         atomic_add(idat.fcy, gy, i);
         atomic_add(idat.fcz, gz, i);
         atomic_add(data[threadIdx.x].fcx, gx, shk);
         atomic_add(data[threadIdx.x].fcy, gy, shk);
         atomic_add(data[threadIdx.x].fcz, gz, shk);
      }
      if CONSTEXPR (do_v) {
         atomic_add(vctlxx, vctlyx, vctlzx, vctlyy, vctlzy, vctlzz, vbuf,
                    offset);
      }
   } // end for (iw)
}


//====================================================================//


template <class Ver, class ETYP, class RADRULE, class EPSRULE, int NOUT>
void echglj_cu3()
{
   const auto& st = *cspatial_unit;
   size_t bufsize = buffer_size();


   real eccut, ecoff;
   real aewald = 0;
   if CONSTEXPR (eq<ETYP, EWALD>()) {
      ecoff = switch_off(switch_ewald);
      // eccut = ecoff; // not used
      PMEUnit pu = epme_unit;
      aewald = pu->aewald;
   } else {
      ecoff = switch_off(switch_charge);
      eccut = switch_cut(switch_charge);
   }
   real f = electric / dielec;


   real evoff = switch_off(switch_vdw);
   real evcut = switch_cut(switch_vdw);
   auto i12 = couple_i12;
   if (vdw_exclude_bond == false)
      i12 = nullptr;


   // Record event `echglj_start` when other kernels on `nonblk` have ended.
   check_rt(cudaEventRecord(echglj_start, nonblk));
   // `echglj_stream` will wait until event `echglj_start` is recorded.
   check_rt(cudaStreamWaitEvent(echglj_stream, echglj_start, 0));


   int nparallel = WARP_SIZE * st.niak;
   int grid1 = (nparallel + CHGLJ4_BDIM - 1) / CHGLJ4_BDIM;
   const auto& attr = get_device_attributes()[idevice];
   int maxthreads = attr.max_threads_per_multiprocessor;
   int mpcount = attr.multiprocessor_count;
   int grid2 = maxthreads / CHGLJ4_BDIM * mpcount;
   int ngrid = std::min(grid1, grid2);


   if CONSTEXPR (eq<ETYP, EWALD>()) {
      if (NOUT == 2 and st.niak > 0) {
         auto ker1 = echglj_cu1<Ver, EWALD, RADRULE, EPSRULE>;
         launch_k1s(nonblk, WARP_SIZE * st.niak, ker1, //
                    bufsize, ec, vir_ec, decx, decy, decz, eccut, ecoff,
                    ebuffer, f, pchg, //
                    ev, vir_ev, devx, devy, devz, atom_rad, atom_eps, evcut,
                    evoff,             //
                    TINKER_IMAGE_ARGS, //
                    st.sorted, st.niak, st.iak, st.lst, n, aewald, i12);
      } else if (NOUT == 1 and st.niak > 0) {
         auto ker1 = echglj_cu4<Ver, EWALD, RADRULE, EPSRULE>;
         ker1<<<ngrid, CHGLJ4_BDIM, 0, echglj_stream>>>(
            bufsize, ec, vir_ec, decx, decy, decz, //
            eccut, ecoff, ebuffer, f, pchg,        //
            atom_rad, atom_eps, evcut, evoff,      //
            TINKER_IMAGE_ARGS,                     //
            st.sorted, st.niak, st.iak, st.lst, n, aewald, i12);
         check_rt(cudaEventRecord(echglj_event, echglj_stream));
      }
      if (ncvexclude > 0) {
         auto ker2 = echglj_cu2<Ver, NON_EWALD, RADRULE, EPSRULE>;
         launch_k1s(nonblk, ncvexclude, ker2, //
                    bufsize, ec, vir_ec, decx, decy, decz, eccut, ecoff,
                    ebuffer, f, pchg, //
                    ev, vir_ev, devx, devy, devz, atom_rad, atom_eps, evcut,
                    evoff,             //
                    TINKER_IMAGE_ARGS, //
                    x, y, z, ncvexclude, cvexclude, cvexclude_scale);
      }
   } else if CONSTEXPR (eq<ETYP, NON_EWALD_TAPER>()) {
      if (NOUT == 2 and st.niak > 0) {
         auto ker1 = echglj_cu1<Ver, NON_EWALD_TAPER, RADRULE, EPSRULE>;
         launch_k1s(nonblk, WARP_SIZE * st.niak, ker1, //
                    bufsize, ec, vir_ec, decx, decy, decz, eccut, ecoff,
                    ebuffer, f, pchg, //
                    ev, vir_ev, devx, devy, devz, atom_rad, atom_eps, evcut,
                    evoff,             //
                    TINKER_IMAGE_ARGS, //
                    st.sorted, st.niak, st.iak, st.lst, n, 0, i12);
      } else if (NOUT == 1 and st.niak > 0) {
         auto ker1 = echglj_cu4<Ver, NON_EWALD_TAPER, RADRULE, EPSRULE>;
         ker1<<<ngrid, CHGLJ4_BDIM, 0, echglj_stream>>>(
            bufsize, ec, vir_ec, decx, decy, decz, //
            eccut, ecoff, ebuffer, f, pchg,        //
            atom_rad, atom_eps, evcut, evoff,      //
            TINKER_IMAGE_ARGS,                     //
            st.sorted, st.niak, st.iak, st.lst, n, 0, i12);
         check_rt(cudaEventRecord(echglj_event, echglj_stream));
      }
      if (ncvexclude > 0) {
         auto ker2 = echglj_cu2<Ver, NON_EWALD_TAPER, RADRULE, EPSRULE>;
         launch_k1s(nonblk, ncvexclude, ker2, //
                    bufsize, ec, vir_ec, decx, decy, decz, eccut, ecoff,
                    ebuffer, f, pchg, //
                    ev, vir_ev, devx, devy, devz, atom_rad, atom_eps, evcut,
                    evoff,             //
                    TINKER_IMAGE_ARGS, //
                    x, y, z, ncvexclude, cvexclude, cvexclude_scale);
      }
   }
}


//====================================================================//


void echglj_rad_arith_eps_geom_nonewald_cu(int vers)
{
   if (use_osrw) {
      constexpr int N = 2;
      if (vers == calc::v0)
         echglj_cu3<calc::V0, NON_EWALD_TAPER, RAD_ARITH, EPS_GEOM, N>();
      else if (vers == calc::v1)
         echglj_cu3<calc::V1, NON_EWALD_TAPER, RAD_ARITH, EPS_GEOM, N>();
      else if (vers == calc::v3)
         echglj_cu3<calc::V3, NON_EWALD_TAPER, RAD_ARITH, EPS_GEOM, N>();
      else if (vers == calc::v4)
         echglj_cu3<calc::V4, NON_EWALD_TAPER, RAD_ARITH, EPS_GEOM, N>();
      else if (vers == calc::v5)
         echglj_cu3<calc::V5, NON_EWALD_TAPER, RAD_ARITH, EPS_GEOM, N>();
      else if (vers == calc::v6)
         echglj_cu3<calc::V6, NON_EWALD_TAPER, RAD_ARITH, EPS_GEOM, N>();
   } else {
      constexpr int N = 1;
      if (vers == calc::v0)
         echglj_cu3<calc::V0, NON_EWALD_TAPER, RAD_ARITH, EPS_GEOM, N>();
      else if (vers == calc::v1)
         echglj_cu3<calc::V1, NON_EWALD_TAPER, RAD_ARITH, EPS_GEOM, N>();
      else if (vers == calc::v3)
         echglj_cu3<calc::V3, NON_EWALD_TAPER, RAD_ARITH, EPS_GEOM, N>();
      else if (vers == calc::v4)
         echglj_cu3<calc::V4, NON_EWALD_TAPER, RAD_ARITH, EPS_GEOM, N>();
      else if (vers == calc::v5)
         echglj_cu3<calc::V5, NON_EWALD_TAPER, RAD_ARITH, EPS_GEOM, N>();
      else if (vers == calc::v6)
         echglj_cu3<calc::V6, NON_EWALD_TAPER, RAD_ARITH, EPS_GEOM, N>();
   }

   elj14_cu(vers);
}


void echglj_rad_arith_eps_geom_ewald_real_cu(int vers)
{
   if (use_osrw) {
      constexpr int N = 2;
      if (vers == calc::v0)
         echglj_cu3<calc::V0, EWALD, RAD_ARITH, EPS_GEOM, N>();
      else if (vers == calc::v1)
         echglj_cu3<calc::V1, EWALD, RAD_ARITH, EPS_GEOM, N>();
      else if (vers == calc::v3)
         echglj_cu3<calc::V3, EWALD, RAD_ARITH, EPS_GEOM, N>();
      else if (vers == calc::v4)
         echglj_cu3<calc::V4, EWALD, RAD_ARITH, EPS_GEOM, N>();
      else if (vers == calc::v5)
         echglj_cu3<calc::V5, EWALD, RAD_ARITH, EPS_GEOM, N>();
      else if (vers == calc::v6)
         echglj_cu3<calc::V6, EWALD, RAD_ARITH, EPS_GEOM, N>();
   } else {
      constexpr int N = 1;
      if (vers == calc::v0)
         echglj_cu3<calc::V0, EWALD, RAD_ARITH, EPS_GEOM, N>();
      else if (vers == calc::v1)
         echglj_cu3<calc::V1, EWALD, RAD_ARITH, EPS_GEOM, N>();
      else if (vers == calc::v3)
         echglj_cu3<calc::V3, EWALD, RAD_ARITH, EPS_GEOM, N>();
      else if (vers == calc::v4)
         echglj_cu3<calc::V4, EWALD, RAD_ARITH, EPS_GEOM, N>();
      else if (vers == calc::v5)
         echglj_cu3<calc::V5, EWALD, RAD_ARITH, EPS_GEOM, N>();
      else if (vers == calc::v6)
         echglj_cu3<calc::V6, EWALD, RAD_ARITH, EPS_GEOM, N>();
   }


   elj14_cu(vers);
}
}
