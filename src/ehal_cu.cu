#include "add.h"
#include "evdw.h"
#include "glob.spatial.h"
#include "image.h"
#include "launch.h"
#include "md.h"
#include "seq_pair_hal.h"
#include "seq_switch.h"
#include "seq_triangle.h"
#include "spatial2.h"
#include "switch.h"
#include "tool/gpu_card.h"


/**
 * Overheads:
 *    - Different vcouple methods.
 *    - PBC type in image().
 *    - Random access to the "i" parameters and gradients.
 *    - (If not hard-coded) ghal, dhal, scexp, scalpha.
 */


/**
 * Kernel ehal_cu1 on GTX 1070 for DHFR2
 * 7 angstroms cutoff and 10 % buffer
 *
 * unsorted (a) | generic image (b) | decouple vlambda (c) | 1e-6 s
 * --------------------------------------------------------------------
 * -            | -                 | -                    | 342
 * +            | -                 | -                    | 348
 * +            | -                 | +                    | 366
 * +            | +                 | -                    | 373
 * +            | +                 | +                    | 382
 * --------------------------------------------------------------------
 * (a) - assuming sorted[i].unsorted == i, no random memory access for the "i"
 *     parameters and gradients; + the original code.
 * (b) - hard-coded orthogonal image routine; + generic image routine.
 * (c) - hard-coded decouple method; + generic vlambda method.
 */


namespace tinker {
#if 1
#   define GHAL    (real)0.12
#   define DHAL    (real)0.07
#   define SCEXP   5
#   define SCALPHA (real)0.7
#elif 0
#   define GHAL    ghal
#   define DHAL    dhal
#   define SCEXP   scexp
#   define SCALPHA scalpha
#endif
template <class Ver>
__global__
void ehal_cu1(count_buffer restrict nebuf, energy_buffer restrict ebuf,
              virial_buffer restrict vbuf, grad_prec* restrict gx,
              grad_prec* restrict gy, grad_prec* restrict gz,
              TINKER_IMAGE_PARAMS, real cut, real off, const real* restrict x,
              const real* restrict y, const real* restrict z, int n,
              const Spatial::SortedAtom* restrict sorted, int nakpl,
              const int* restrict iakpl, int niak, const int* restrict iak,
              const int* restrict lst, int nexclude,
              const int (*restrict exclude)[2],
              const real* restrict exclude_scale,
              const unsigned int* restrict info, int njvdw, real vlam,
              evdw_t vcouple, const real* restrict radmin,
              const real* restrict epsilon, const int* restrict jvdw,
              const int* restrict mut)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);


   using ebuf_prec = energy_buffer_traits::type;
   using vbuf_prec = virial_buffer_traits::type;
   ebuf_prec etl;
   vbuf_prec vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz;
   if CONSTEXPR (do_e) {
      etl = 0;
   }
   if CONSTEXPR (do_v) {
      vtlxx = 0;
      vtlyx = 0;
      vtlzx = 0;
      vtlyy = 0;
      vtlzy = 0;
      vtlzz = 0;
   }
   int ctl;
   if CONSTEXPR (do_a) {
      ctl = 0;
   }


   //* /
   // exclude
   for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
      int i = exclude[ii][0];
      int k = exclude[ii][1];
      real scale = exclude_scale[ii];


      int ijvdw = jvdw[i];
      int imut = mut[i];
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];


      int kjvdw = jvdw[k];
      int kmut = mut[k];
      real xr = xi - x[k];
      real yr = yi - y[k];
      real zr = zi - z[k];


      real r2 = image2(xr, yr, zr);
      real r = REAL_SQRT(r2);
      real invr = REAL_RECIP(r);
      real rv = radmin[ijvdw * njvdw + kjvdw];
      real eps = epsilon[ijvdw * njvdw + kjvdw];
      real vlambda = 1;
      if (vcouple == evdw_t::decouple) {
         vlambda = (imut == kmut ? 1 : vlam);
      } else if (vcouple == evdw_t::annihilate) {
         vlambda = (imut || kmut ? vlam : 1);
      }
      real e, de;
      pair_hal_v2<do_g, 0>(r, scale, rv, eps, cut, off, vlambda, GHAL, DHAL,
                           SCEXP, SCALPHA, e, de);


      if CONSTEXPR (do_e) {
         etl += cvt_to<ebuf_prec>(e);
         if CONSTEXPR (do_a) {
            if (scale != 0 and e != 0) {
               ctl += 1;
            }
         }
      }
      if CONSTEXPR (do_g) {
         real dedx, dedy, dedz;
         de *= invr;
         dedx = de * xr;
         dedy = de * yr;
         dedz = de * zr;
         atomic_add(dedx, gx, i);
         atomic_add(dedy, gy, i);
         atomic_add(dedz, gz, i);
         atomic_add(-dedx, gx, k);
         atomic_add(-dedy, gy, k);
         atomic_add(-dedz, gz, k);
         if CONSTEXPR (do_v) {
            vtlxx += cvt_to<vbuf_prec>(xr * dedx);
            vtlyx += cvt_to<vbuf_prec>(yr * dedx);
            vtlzx += cvt_to<vbuf_prec>(zr * dedx);
            vtlyy += cvt_to<vbuf_prec>(yr * dedy);
            vtlzy += cvt_to<vbuf_prec>(zr * dedy);
            vtlzz += cvt_to<vbuf_prec>(zr * dedz);
         }
      }
   }
   // */


   //* /
   // block pairs that have scale factors
   for (int iw = iwarp; iw < nakpl; iw += nwarp) {
      real fix, fiy, fiz;
      real fkx, fky, fkz;
      if CONSTEXPR (do_g) {
         fix = 0;
         fiy = 0;
         fiz = 0;
         fkx = 0;
         fky = 0;
         fkz = 0;
      }


      int tri, tx, ty;
      tri = iakpl[iw];
      tri_to_xy(tri, tx, ty);


      int shiid = ty * WARP_SIZE + ilane;
      int shatomi = min(shiid, n - 1);
      real shxi = sorted[shatomi].x;
      real shyi = sorted[shatomi].y;
      real shzi = sorted[shatomi].z;
      int shi = sorted[shatomi].unsorted;
      int shijvdw = jvdw[shi];
      int shimut = mut[shi];


      int kid = tx * WARP_SIZE + ilane;
      int atomk = min(kid, n - 1);
      real xk = sorted[atomk].x;
      real yk = sorted[atomk].y;
      real zk = sorted[atomk].z;
      int k = sorted[atomk].unsorted;
      int kjvdw = jvdw[k];
      int kmut = mut[k];


      int bit0 = info[iw * WARP_SIZE + ilane];
      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int srcmask = 1 << srclane;
         int bit = bit0 & srcmask;
         int ijvdw = shijvdw;
         int imut = shimut;
         int iid = shiid;
         real xr = shxi - xk;
         real yr = shyi - yk;
         real zr = shzi - zk;


         bool incl = iid < kid and kid < n and bit == 0;
         real r2 = image2(xr, yr, zr);
         real r = REAL_SQRT(r2);
         real invr = REAL_RECIP(r);
         real rv = radmin[ijvdw * njvdw + kjvdw];
         real eps = epsilon[ijvdw * njvdw + kjvdw];
         real vlambda = 1;
         if (vcouple == evdw_t::decouple) {
            vlambda = (imut == kmut ? 1 : vlam);
         } else if (vcouple == evdw_t::annihilate) {
            vlambda = (imut || kmut ? vlam : 1);
         }
         real e, de;
         pair_hal_v2<do_g, 1>(r, 1, rv, eps, cut, off, vlambda, GHAL, DHAL,
                              SCEXP, SCALPHA, e, de);


         if CONSTEXPR (do_e) {
            etl += incl ? cvt_to<ebuf_prec>(e) : 0;
            if CONSTEXPR (do_a) {
               if (incl and e != 0) {
                  ctl += 1;
               }
            }
         }
         if CONSTEXPR (do_g) {
            real dedx, dedy, dedz;
            de = incl ? de * invr : 0;
            dedx = de * xr;
            dedy = de * yr;
            dedz = de * zr;
            fix += dedx;
            fiy += dedy;
            fiz += dedz;
            fkx -= dedx;
            fky -= dedy;
            fkz -= dedz;
            if CONSTEXPR (do_v) {
               vtlxx += cvt_to<vbuf_prec>(xr * dedx);
               vtlyx += cvt_to<vbuf_prec>(yr * dedx);
               vtlzx += cvt_to<vbuf_prec>(zr * dedx);
               vtlyy += cvt_to<vbuf_prec>(yr * dedy);
               vtlzy += cvt_to<vbuf_prec>(zr * dedy);
               vtlzz += cvt_to<vbuf_prec>(zr * dedz);
            }
         }


         shijvdw = __shfl_sync(ALL_LANES, shijvdw, ilane + 1);
         shimut = __shfl_sync(ALL_LANES, shimut, ilane + 1);
         shiid = __shfl_sync(ALL_LANES, shiid, ilane + 1);
         shxi = __shfl_sync(ALL_LANES, shxi, ilane + 1);
         shyi = __shfl_sync(ALL_LANES, shyi, ilane + 1);
         shzi = __shfl_sync(ALL_LANES, shzi, ilane + 1);
         fix = __shfl_sync(ALL_LANES, fix, ilane + 1);
         fiy = __shfl_sync(ALL_LANES, fiy, ilane + 1);
         fiz = __shfl_sync(ALL_LANES, fiz, ilane + 1);
      }


      if CONSTEXPR (do_g) {
         atomic_add(fix, gx, shi);
         atomic_add(fiy, gy, shi);
         atomic_add(fiz, gz, shi);
         atomic_add(fkx, gx, k);
         atomic_add(fky, gy, k);
         atomic_add(fkz, gz, k);
      }
   } // end loop block pairs
   // */


   //* /
   // block-atoms
   for (int iw = iwarp; iw < niak; iw += nwarp) {
      real fix, fiy, fiz;
      real fkx, fky, fkz;
      if CONSTEXPR (do_g) {
         fix = 0;
         fiy = 0;
         fiz = 0;
         fkx = 0;
         fky = 0;
         fkz = 0;
      }


      int ty = iak[iw];
      int shatomi = ty * WARP_SIZE + ilane;
      real shxi = sorted[shatomi].x;
      real shyi = sorted[shatomi].y;
      real shzi = sorted[shatomi].z;
      int shi = sorted[shatomi].unsorted;
      int shijvdw = jvdw[shi];
      int shimut = mut[shi];


      int atomk = lst[iw * WARP_SIZE + ilane];
      real xk = sorted[atomk].x;
      real yk = sorted[atomk].y;
      real zk = sorted[atomk].z;
      int k = sorted[atomk].unsorted;
      int kjvdw = jvdw[k];
      int kmut = mut[k];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int ijvdw = shijvdw;
         int imut = shimut;
         real xr = shxi - xk;
         real yr = shyi - yk;
         real zr = shzi - zk;


         bool incl = atomk > 0;
         real r2 = image2(xr, yr, zr);
         real r = REAL_SQRT(r2);
         real invr = REAL_RECIP(r);
         real rv = radmin[ijvdw * njvdw + kjvdw];
         real eps = epsilon[ijvdw * njvdw + kjvdw];
         real vlambda = 1;
         if (vcouple == evdw_t::decouple) {
            vlambda = (imut == kmut ? 1 : vlam);
         } else if (vcouple == evdw_t::annihilate) {
            vlambda = (imut || kmut ? vlam : 1);
         }
         real e, de;
         pair_hal_v2<do_g, 1>(r, 1, rv, eps, cut, off, vlambda, GHAL, DHAL,
                              SCEXP, SCALPHA, e, de);


         if CONSTEXPR (do_e) {
            etl += incl ? cvt_to<ebuf_prec>(e) : 0;
            if CONSTEXPR (do_a) {
               if (incl and e != 0) {
                  ctl += 1;
               }
            }
         }
         if CONSTEXPR (do_g) {
            real dedx, dedy, dedz;
            de = incl ? de * invr : 0;
            dedx = de * xr;
            dedy = de * yr;
            dedz = de * zr;
            fix += dedx;
            fiy += dedy;
            fiz += dedz;
            fkx -= dedx;
            fky -= dedy;
            fkz -= dedz;
            if CONSTEXPR (do_v) {
               vtlxx += cvt_to<vbuf_prec>(xr * dedx);
               vtlyx += cvt_to<vbuf_prec>(yr * dedx);
               vtlzx += cvt_to<vbuf_prec>(zr * dedx);
               vtlyy += cvt_to<vbuf_prec>(yr * dedy);
               vtlzy += cvt_to<vbuf_prec>(zr * dedy);
               vtlzz += cvt_to<vbuf_prec>(zr * dedz);
            }
         }


         shijvdw = __shfl_sync(ALL_LANES, shijvdw, ilane + 1);
         shimut = __shfl_sync(ALL_LANES, shimut, ilane + 1);
         shxi = __shfl_sync(ALL_LANES, shxi, ilane + 1);
         shyi = __shfl_sync(ALL_LANES, shyi, ilane + 1);
         shzi = __shfl_sync(ALL_LANES, shzi, ilane + 1);
         fix = __shfl_sync(ALL_LANES, fix, ilane + 1);
         fiy = __shfl_sync(ALL_LANES, fiy, ilane + 1);
         fiz = __shfl_sync(ALL_LANES, fiz, ilane + 1);
      }


      if CONSTEXPR (do_g) {
         atomic_add(fix, gx, shi);
         atomic_add(fiy, gy, shi);
         atomic_add(fiz, gz, shi);
         atomic_add(fkx, gx, k);
         atomic_add(fky, gy, k);
         atomic_add(fkz, gz, k);
      }
   } // end loop block-atoms
   // */


   if CONSTEXPR (do_a) {
      atomic_add(ctl, nebuf, ithread);
   }
   if CONSTEXPR (do_e) {
      atomic_add(etl, ebuf, ithread);
   }
   if CONSTEXPR (do_v) {
      atomic_add(vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz, vbuf, ithread);
   }
}


template <class Ver>
void ehal_cu3()
{
   constexpr bool do_g = Ver::g;


   const auto& st = *vspatial_v2_unit;
   const real cut = switch_cut(switch_vdw);
   const real off = switch_off(switch_vdw);


   if CONSTEXPR (do_g)
      darray::zero(PROCEED_NEW_Q, n, gxred, gyred, gzred);


   int ngrid = get_grid_size(BLOCK_DIM);
   auto ker1 = ehal_cu1<Ver>;
   ker1<<<ngrid, BLOCK_DIM, 0, nonblk>>>(
      nev, ev, vir_ev, gxred, gyred, gzred, TINKER_IMAGE_ARGS, cut, off, st.x,
      st.y, st.z, st.n, st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst,
      nvexclude, vexclude, vexclude_scale, st.si1.bit0, njvdw, vlam, vcouple,
      radmin, epsilon, jvdw, mut);


   if CONSTEXPR (do_g) {
      ehal_resolve_gradient();
   }
}


void ehal_cu(int vers)
{
   if (vers == calc::v0)
      ehal_cu3<calc::V0>();
   else if (vers == calc::v1)
      ehal_cu3<calc::V1>();
   else if (vers == calc::v3)
      ehal_cu3<calc::V3>();
   else if (vers == calc::v4)
      ehal_cu3<calc::V4>();
   else if (vers == calc::v5)
      ehal_cu3<calc::V5>();
   else if (vers == calc::v6)
      ehal_cu3<calc::V6>();
}
}
