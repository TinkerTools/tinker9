#include "add.h"
#include "echarge.h"
#include "image.h"
#include "launch.h"
#include "mdegv.h"
#include "mdpq.h"
#include "named_struct.h"
#include "pmestuf.h"
#include "seq_bsplgen.h"
#include "seq_pair_charge.h"
#include "spatial.h"
#include "switch.h"


namespace tinker {
#define ECHARGE_ARGS                                                           \
   size_t bufsize, count_buffer restrict nec, energy_buffer restrict ec,       \
      virial_buffer restrict vir_ec, grad_prec *restrict gx,                   \
      grad_prec *restrict gy, grad_prec *restrict gz, TINKER_IMAGE_PARAMS,     \
      real cut, real off, real ebuffer, real f, const real *restrict pchg


template <class Ver, class ETYP>
__global__
void echarge_cu1(ECHARGE_ARGS, const Spatial::SortedAtom* restrict sorted,
                 int niak, const int* restrict iak, const int* restrict lst,
                 int n, real aewald)
{
   static_assert(eq<ETYP, EWALD>() || eq<ETYP, NON_EWALD_TAPER>(), "");
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   const int iwarp = (threadIdx.x + blockIdx.x * blockDim.x) / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);
   const int offset = (threadIdx.x + blockIdx.x * blockDim.x) & (bufsize - 1);


   struct Data
   {
      real x, y, z, chg;
      real frcx, frcy, frcz;
      real padding_;
   };
   __shared__ Data data[BLOCK_DIM];


   const real off2 = off * off;
   for (int iw = iwarp; iw < niak; iw += nwarp) {
      int ctl;
      if CONSTEXPR (do_a) {
         ctl = 0;
      }
      real etl;
      if CONSTEXPR (do_e) {
         etl = 0;
      }
      real vtlxx, vtlxy, vtlxz, vtlyy, vtlyz, vtlzz;
      if CONSTEXPR (do_v) {
         vtlxx = 0;
         vtlxy = 0;
         vtlxz = 0;
         vtlyy = 0;
         vtlyz = 0;
         vtlzz = 0;
      }


      Data idat;
      if CONSTEXPR (do_g) {
         idat.frcx = 0;
         idat.frcy = 0;
         idat.frcz = 0;
      }
      int atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
      idat.x = sorted[atomi].x;
      idat.y = sorted[atomi].y;
      idat.z = sorted[atomi].z;
      int i = sorted[atomi].unsorted;
      idat.chg = pchg[i];


      if CONSTEXPR (do_g) {
         data[threadIdx.x].frcx = 0;
         data[threadIdx.x].frcy = 0;
         data[threadIdx.x].frcz = 0;
      }
      int shatomk = lst[iw * WARP_SIZE + ilane];
      data[threadIdx.x].x = sorted[shatomk].x;
      data[threadIdx.x].y = sorted[shatomk].y;
      data[threadIdx.x].z = sorted[shatomk].z;
      int shk = sorted[shatomk].unsorted;
      data[threadIdx.x].chg = pchg[shk];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
         real xr = idat.x - data[klane].x;
         real yr = idat.y - data[klane].y;
         real zr = idat.z - data[klane].z;


         real r2 = image2(xr, yr, zr);
         if (atomi < atomk && r2 <= off2) {
            real r = REAL_SQRT(r2);


            MAYBE_UNUSED real grdx, grdy, grdz;
            if CONSTEXPR (do_g) {
               grdx = 0;
               grdy = 0;
               grdz = 0;
            }


            pair_charge<Ver, ETYP>(r, xr, yr, zr, 1, idat.chg, data[klane].chg,
                                   ebuffer, f, aewald, cut, off, //
                                   grdx, grdy, grdz, ctl, etl, vtlxx, vtlxy,
                                   vtlxz, vtlyy, vtlyz, vtlzz);


            if CONSTEXPR (do_g) {
               idat.frcx += grdx;
               idat.frcy += grdy;
               idat.frcz += grdz;
               data[klane].frcx -= grdx;
               data[klane].frcy -= grdy;
               data[klane].frcz -= grdz;
            }
         } // end if (include)
      }


      if CONSTEXPR (do_a)
         atomic_add(ctl, nec, offset);
      if CONSTEXPR (do_e)
         atomic_add(etl, ec, offset);
      if CONSTEXPR (do_g) {
         atomic_add(idat.frcx, &gx[i]);
         atomic_add(idat.frcy, &gy[i]);
         atomic_add(idat.frcz, &gz[i]);
         atomic_add(data[threadIdx.x].frcx, &gx[shk]);
         atomic_add(data[threadIdx.x].frcy, &gy[shk]);
         atomic_add(data[threadIdx.x].frcz, &gz[shk]);
      }
      if CONSTEXPR (do_v)
         atomic_add(vtlxx, vtlxy, vtlxz, vtlyy, vtlyz, vtlzz, vir_ec, offset);
   } // end for (iw)
}


template <class Ver, class ETYP>
__global__
void echarge_cu2(ECHARGE_ARGS, const real* restrict x, const real* restrict y,
                 const real* restrict z, int ncexclude,
                 const int (*restrict cexclude)[2],
                 const real* restrict cexclude_scale)
{
   static_assert(eq<ETYP, NON_EWALD>() || eq<ETYP, NON_EWALD_TAPER>(), "");
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < ncexclude;
        ii += blockDim.x * gridDim.x) {
      int offset = ii & (bufsize - 1);


      int i = cexclude[ii][0];
      int k = cexclude[ii][1];
      real cscale = cexclude_scale[ii];


      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real ci = pchg[i];


      real xr = xi - x[k];
      real yr = yi - y[k];
      real zr = zi - z[k];
      real ck = pchg[k];


      real r2 = image2(xr, yr, zr);
      real off2 = off * off;
      if (r2 <= off2) {
         MAYBE_UNUSED int ctl;
         MAYBE_UNUSED real e, grdx, grdy, grdz, vxx, vxy, vxz, vyy, vyz, vzz;
         if CONSTEXPR (do_a) {
            ctl = 0;
         }
         if CONSTEXPR (do_e) {
            e = 0;
         }
         if CONSTEXPR (do_g) {
            grdx = 0;
            grdy = 0;
            grdz = 0;
         }
         if CONSTEXPR (do_v) {
            vxx = 0;
            vxy = 0;
            vxz = 0;
            vyy = 0;
            vyz = 0;
            vzz = 0;
         }


         real r = REAL_SQRT(r2);
         pair_charge<Ver, ETYP>(
            r, xr, yr, zr, cscale, ci, ck, ebuffer, f, 0, cut, off, //
            grdx, grdy, grdz, ctl, e, vxx, vxy, vxz, vyy, vyz, vzz);
         if (e != 0) {
            if CONSTEXPR (do_a)
               atomic_add(ctl, nec, offset);
            if CONSTEXPR (do_e)
               atomic_add(e, ec, offset);
         }
         if CONSTEXPR (do_g) {
            atomic_add(grdx, gx, i);
            atomic_add(grdy, gy, i);
            atomic_add(grdz, gz, i);
            atomic_add(-grdx, gx, k);
            atomic_add(-grdy, gy, k);
            atomic_add(-grdz, gz, k);
         }
         if CONSTEXPR (do_v)
            atomic_add(vxx, vxy, vxz, vyy, vyz, vzz, vir_ec, offset);
      } // end if (include)
   }
}


template <class Ver, class ETYP>
void echarge_cu()
{
   const auto& st = *cspatial_unit;
   real cut, off;
   if CONSTEXPR (eq<ETYP, EWALD>()) {
      off = switch_off(switch_ewald);
      // cut = off; // not used
   } else {
      off = switch_off(switch_charge);
      cut = switch_cut(switch_charge);
   }
   auto bufsize = buffer_size();


   const real f = electric / dielec;
   real aewald = 0;
   if CONSTEXPR (eq<ETYP, EWALD>()) {
      PMEUnit pu = epme_unit;
      aewald = pu->aewald;
   }


   if CONSTEXPR (eq<ETYP, EWALD>()) {
      if (st.niak > 0) {
         auto ker1 = echarge_cu1<Ver, EWALD>;
         launch_k1s(nonblk, WARP_SIZE * st.niak, ker1, //
                    bufsize, nec, ec, vir_ec, decx, decy, decz,
                    TINKER_IMAGE_ARGS, 0, off, ebuffer, f, pchg, //
                    st.sorted, st.niak, st.iak, st.lst, n, aewald);
      }
      if (ncexclude > 0) {
         auto ker2 = echarge_cu2<Ver, NON_EWALD>;
         launch_k1s(nonblk, ncexclude, ker2, //
                    bufsize, nec, ec, vir_ec, decx, decy, decz,
                    TINKER_IMAGE_ARGS, 0, off, ebuffer, f, pchg, //
                    x, y, z, ncexclude, cexclude, cexclude_scale);
      }
   } else if CONSTEXPR (eq<ETYP, NON_EWALD_TAPER>()) {
      if (st.niak > 0) {
         auto ker1 = echarge_cu1<Ver, NON_EWALD_TAPER>;
         launch_k1s(nonblk, WARP_SIZE * st.niak, ker1, //
                    bufsize, nec, ec, vir_ec, decx, decy, decz,
                    TINKER_IMAGE_ARGS, cut, off, ebuffer, f, pchg, //
                    st.sorted, st.niak, st.iak, st.lst, n, 0);
      }
      if (ncexclude > 0) {
         auto ker2 = echarge_cu2<Ver, NON_EWALD_TAPER>;
         launch_k1s(nonblk, ncexclude, ker2, //
                    bufsize, nec, ec, vir_ec, decx, decy, decz,
                    TINKER_IMAGE_ARGS, cut, off, ebuffer, f, pchg, //
                    x, y, z, ncexclude, cexclude, cexclude_scale);
      }
   }
}


void echarge_nonewald_cu(int vers)
{
   if (vers == calc::v0)
      echarge_cu<calc::V0, NON_EWALD_TAPER>();
   else if (vers == calc::v1)
      echarge_cu<calc::V1, NON_EWALD_TAPER>();
   else if (vers == calc::v3)
      echarge_cu<calc::V3, NON_EWALD_TAPER>();
   else if (vers == calc::v4)
      echarge_cu<calc::V4, NON_EWALD_TAPER>();
   else if (vers == calc::v5)
      echarge_cu<calc::V5, NON_EWALD_TAPER>();
   else if (vers == calc::v6)
      echarge_cu<calc::V6, NON_EWALD_TAPER>();
}


void echarge_ewald_real_cu(int vers)
{
   if (vers == calc::v0)
      echarge_cu<calc::V0, EWALD>();
   else if (vers == calc::v1)
      echarge_cu<calc::V1, EWALD>();
   else if (vers == calc::v3)
      echarge_cu<calc::V3, EWALD>();
   else if (vers == calc::v4)
      echarge_cu<calc::V4, EWALD>();
   else if (vers == calc::v5)
      echarge_cu<calc::V5, EWALD>();
   else if (vers == calc::v6)
      echarge_cu<calc::V6, EWALD>();
}


//====================================================================//


template <class Ver, int bsorder>
__global__
void echarge_cu3(size_t bufsize, count_buffer restrict nec,
                 energy_buffer restrict ec, const real* restrict pchg, real f,
                 real aewald, int n, int nfft1, int nfft2, int nfft3,
                 const real* restrict x, const real* restrict y,
                 const real* restrict z, const real* restrict qgrid, real3 reca,
                 real3 recb, real3 recc, grad_prec* restrict gx,
                 grad_prec* restrict gy, grad_prec* restrict gz)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;


   real thetai1[4 * 5];
   real thetai2[4 * 5];
   real thetai3[4 * 5];
   __shared__ real sharedarray[5 * 5 * PME_BLOCKDIM];
   real* restrict array = &sharedarray[5 * 5 * threadIdx.x];


   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < n;
        ii += blockDim.x * gridDim.x) {
      real chgi = pchg[ii];
      if (chgi == 0)
         continue;


      // self energy, tinfoil
      if CONSTEXPR (do_e) {
         int offset = ii & (bufsize - 1);
         real fs = -f * aewald * REAL_RECIP(sqrtpi);
         real e = fs * chgi * chgi;
         atomic_add(e, ec, offset);
         if (do_a) {
            atomic_add(1, nec, offset);
         }
      }


      // recip gradient
      if CONSTEXPR (do_g) {
         real xi = x[ii];
         real yi = y[ii];
         real zi = z[ii];


         real w1 = xi * reca.x + yi * reca.y + zi * reca.z;
         w1 = w1 + 0.5f - REAL_FLOOR(w1 + 0.5f);
         real fr1 = nfft1 * w1;
         int igrid1 = REAL_FLOOR(fr1);
         w1 = fr1 - igrid1;


         real w2 = xi * recb.x + yi * recb.y + zi * recb.z;
         w2 = w2 + 0.5f - REAL_FLOOR(w2 + 0.5f);
         real fr2 = nfft2 * w2;
         int igrid2 = REAL_FLOOR(fr2);
         w2 = fr2 - igrid2;


         real w3 = xi * recc.x + yi * recc.y + zi * recc.z;
         w3 = w3 + 0.5f - REAL_FLOOR(w3 + 0.5f);
         real fr3 = nfft3 * w3;
         int igrid3 = REAL_FLOOR(fr3);
         w3 = fr3 - igrid3;


         igrid1 = igrid1 - bsorder + 1;
         igrid2 = igrid2 - bsorder + 1;
         igrid3 = igrid3 - bsorder + 1;
         igrid1 += (igrid1 < 0 ? nfft1 : 0);
         igrid2 += (igrid2 < 0 ? nfft2 : 0);
         igrid3 += (igrid3 < 0 ? nfft3 : 0);


         bsplgen<2, bsorder>(w1, thetai1, array);
         bsplgen<2, bsorder>(w2, thetai2, array);
         bsplgen<2, bsorder>(w3, thetai3, array);


         real fi = f * chgi;
         real de1 = 0, de2 = 0, de3 = 0;
         for (int iz = 0; iz < bsorder; ++iz) {
            int zbase = igrid3 + iz;
            zbase -= (zbase >= nfft3 ? nfft3 : 0);
            zbase *= (nfft1 * nfft2);
            real t3 = thetai3[4 * iz];
            real dt3 = nfft3 * thetai3[1 + 4 * iz];
            for (int iy = 0; iy < bsorder; ++iy) {
               int ybase = igrid2 + iy;
               ybase -= (ybase >= nfft2 ? nfft2 : 0);
               ybase *= nfft1;
               real t2 = thetai2[4 * iy];
               real dt2 = nfft2 * thetai2[1 + 4 * iy];
               for (int ix = 0; ix < bsorder; ++ix) {
                  int xbase = igrid1 + ix;
                  xbase -= (xbase >= nfft1 ? nfft1 : 0);
                  int index = xbase + ybase + zbase;
                  real term = qgrid[2 * index];
                  real t1 = thetai1[4 * ix];
                  real dt1 = nfft1 * thetai1[1 + 4 * ix];
                  de1 += term * dt1 * t2 * t3;
                  de2 += term * dt2 * t1 * t3;
                  de3 += term * dt3 * t1 * t2;
               }
            }
         } // end for (iz)


         real frcx = fi * (reca.x * de1 + recb.x * de2 + recc.x * de3);
         real frcy = fi * (reca.y * de1 + recb.y * de2 + recc.y * de3);
         real frcz = fi * (reca.z * de1 + recb.z * de2 + recc.z * de3);
         atomic_add(frcx, gx, ii);
         atomic_add(frcy, gy, ii);
         atomic_add(frcz, gz, ii);
      }
   }
}


template <class Ver>
void echarge_fphi_self_cu()
{
   auto bufsize = buffer_size();
   real f = electric / dielec;
   const auto& st = *epme_unit;
   real aewald = st.aewald;
   int nfft1 = st.nfft1;
   int nfft2 = st.nfft2;
   int nfft3 = st.nfft3;


   auto ker = echarge_cu3<Ver, 5>;
   launch_k2s(nonblk, PME_BLOCKDIM, n, ker,         //
              bufsize, nec, ec, pchg, f, aewald, n, //
              nfft1, nfft2, nfft3, x, y, z, st.qgrid, recipa, recipb, recipc,
              decx, decy, decz);
}


void echarge_ewald_fphi_self_cu(int vers)
{
   if (vers == calc::v0)
      echarge_fphi_self_cu<calc::V0>();
   else if (vers == calc::v1)
      echarge_fphi_self_cu<calc::V1>();
   else if (vers == calc::v3)
      echarge_fphi_self_cu<calc::V3>();
   else if (vers == calc::v4)
      echarge_fphi_self_cu<calc::V4>();
   else if (vers == calc::v5)
      echarge_fphi_self_cu<calc::V5>();
   else if (vers == calc::v6)
      echarge_fphi_self_cu<calc::V6>();
}
}
