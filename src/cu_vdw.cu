#include "add.cuh"
#include "e_vdw.h"
#include "launch.cuh"
#include "md.h"
#include "seq_image.h"
#include "seq_pair_hal.h"
#include "seq_switch.h"
#include "spatial.h"


/**
 * Overheads:
 *    - Different vcouple methods.
 *    - PBC type in image().
 *    - Random access to the "i" parameters and gradients.
 *    - (If not hard-coded) ghal, dhal, scexp, scalpha.
 */


TINKER_NAMESPACE_BEGIN
#define HAL_ARGS                                                               \
   size_t bufsize, count_buffer restrict nev, energy_buffer restrict ev,       \
      virial_buffer restrict vir_ev, real *restrict gxred,                     \
      real *restrict gyred, real *restrict gzred, TINKER_IMAGE_PARAMS,         \
      int njvdw, const int *restrict jvdw, const real *restrict radmin,        \
      const real *restrict epsilon, const real *vlam, evdw_t vcouple,          \
      real cut, real off
#if 1
#   define GHAL (real)0.12
#   define DHAL (real)0.07
#   define SCEXP 5
#   define SCALPHA (real)0.7
#elif 0
#   define GHAL ghal
#   define DHAL dhal
#   define SCEXP scexp
#   define SCALPHA scalpha
#endif
template <int USE>
__launch_bounds__(BLOCK_DIM) __global__
void evdw_hal_cu1(HAL_ARGS, int n, const Spatial::SortedAtom* restrict sorted,
                  int niak, const int* restrict iak, const int* restrict lst)
{
   constexpr int do_e = USE & calc::energy;
   constexpr int do_a = USE & calc::analyz;
   constexpr int do_g = USE & calc::grad;
   constexpr int do_v = USE & calc::virial;


   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);
   const int offset = ithread & (bufsize - 1);


   // thread local variables
   MAYBE_UNUSED int ctl;
   MAYBE_UNUSED real etl;
   MAYBE_UNUSED real gxi, gyi, gzi, gxk, gyk, gzk;
   MAYBE_UNUSED real vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz;


   const real cut2 = cut * cut;
   const real off2 = off * off;
   for (int iw = iwarp; iw < niak; iw += nwarp) {
      if CONSTEXPR (do_a)
         ctl = 0;
      if CONSTEXPR (do_e)
         etl = 0;
      if CONSTEXPR (do_g) {
         gxi = 0;
         gyi = 0;
         gzi = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
      }
      if CONSTEXPR (do_v) {
         vtlxx = 0;
         vtlyx = 0;
         vtlzx = 0;
         vtlyy = 0;
         vtlzy = 0;
         vtlzz = 0;
      }


      int atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
      real xi = sorted[atomi].x;
      real yi = sorted[atomi].y;
      real zi = sorted[atomi].z;
      int i = sorted[atomi].unsorted;
      int it = jvdw[i];
      real lam1 = vlam[i];


      int shatomk = lst[iw * WARP_SIZE + ilane];
      real shx = sorted[shatomk].x;
      real shy = sorted[shatomk].y;
      real shz = sorted[shatomk].z;
      int shk = sorted[shatomk].unsorted;
      int shkt = jvdw[shk];
      real shlam = vlam[shk];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
         real xr = xi - __shfl_sync(ALL_LANES, shx, srclane);
         real yr = yi - __shfl_sync(ALL_LANES, shy, srclane);
         real zr = zi - __shfl_sync(ALL_LANES, shz, srclane);
         int kt = __shfl_sync(ALL_LANES, shkt, srclane);
         real vlambda = __shfl_sync(ALL_LANES, shlam, srclane);


         MAYBE_UNUSED real dedx = 0, dedy = 0, dedz = 0;
         real rik2 = image2(xr, yr, zr);


         if (atomi < atomk && rik2 <= off2) {
            real rik = REAL_SQRT(rik2);
            real rv = radmin[it * njvdw + kt];
            real eps = epsilon[it * njvdw + kt];
            if (vcouple == evdw_t::decouple) {
               vlambda = (lam1 == vlambda ? 1 : REAL_MIN(lam1, vlambda));
            } else if (vcouple == evdw_t::annihilate) {
               vlambda = REAL_MIN(lam1, vlambda);
            }


            MAYBE_UNUSED real e, de;
            pair_hal<do_g>(rik, rv, eps, 1, vlambda, GHAL, DHAL, SCEXP, SCALPHA,
                           e, de);
            if (rik2 > cut2) {
               real taper, dtaper;
               switch_taper5<do_g>(rik, cut, off, taper, dtaper);
               if CONSTEXPR (do_g)
                  de = e * dtaper + de * taper;
               if CONSTEXPR (do_e)
                  e = e * taper;
            }


            if CONSTEXPR (do_a)
               ctl += 1;
            if CONSTEXPR (do_e)
               etl += e;
            if CONSTEXPR (do_g) {
               de *= REAL_RECIP(rik);
               dedx = de * xr;
               dedy = de * yr;
               dedz = de * zr;
               if CONSTEXPR (do_v) {
                  vtlxx += xr * dedx;
                  vtlyx += yr * dedx;
                  vtlzx += zr * dedx;
                  vtlyy += yr * dedy;
                  vtlzy += zr * dedy;
                  vtlzz += zr * dedz;
               }
            }
         } // end if (include)


         if CONSTEXPR (do_g) {
            int dstlane = (ilane + WARP_SIZE - j) & (WARP_SIZE - 1);
            gxi += dedx;
            gyi += dedy;
            gzi += dedz;
            gxk -= __shfl_sync(ALL_LANES, dedx, dstlane);
            gyk -= __shfl_sync(ALL_LANES, dedy, dstlane);
            gzk -= __shfl_sync(ALL_LANES, dedz, dstlane);
         }
      }


      if CONSTEXPR (do_a)
         atomic_add(ctl, nev, offset);
      if CONSTEXPR (do_e)
         atomic_add(etl, ev, offset);
      if CONSTEXPR (do_g) {
         atomic_add(gxi, gxred, i);
         atomic_add(gyi, gyred, i);
         atomic_add(gzi, gzred, i);
         atomic_add(gxk, gxred, shk);
         atomic_add(gyk, gyred, shk);
         atomic_add(gzk, gzred, shk);
      }
      if CONSTEXPR (do_v)
         atomic_add(vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz, vir_ev, offset);
   } // end for (iw)
}


template <int USE>
__global__
void evdw_hal_cu2(HAL_ARGS, const real* restrict xred,
                  const real* restrict yred, const real* restrict zred,
                  int nvexclude_, int (*restrict vexclude_)[2],
                  real* restrict vexclude_scale_)
{
   constexpr int do_e = USE & calc::energy;
   constexpr int do_a = USE & calc::analyz;
   constexpr int do_g = USE & calc::grad;
   constexpr int do_v = USE & calc::virial;


   const real cut2 = cut * cut;
   const real off2 = off * off;
   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < nvexclude_;
        ii += blockDim.x * gridDim.x) {
      int offset = ii & (bufsize - 1);


      int i = vexclude_[ii][0];
      int k = vexclude_[ii][1];
      real vscale = vexclude_scale_[ii];


      int it = jvdw[i];
      real xi = xred[i];
      real yi = yred[i];
      real zi = zred[i];
      real lam1 = vlam[i];


      int kt = jvdw[k];
      real xr = xi - xred[k];
      real yr = yi - yred[k];
      real zr = zi - zred[k];
      real vlambda = vlam[k];


      real rik2 = image2(xr, yr, zr);
      if (rik2 <= off2) {
         real rik = REAL_SQRT(rik2);
         real rv = radmin[it * njvdw + kt];
         real eps = epsilon[it * njvdw + kt];
         if (vcouple == evdw_t::decouple) {
            vlambda = (lam1 == vlambda ? 1 : REAL_MIN(lam1, vlambda));
         } else if (vcouple == evdw_t::annihilate) {
            vlambda = REAL_MIN(lam1, vlambda);
         }


         MAYBE_UNUSED real e, de;
         pair_hal<do_g>(rik, rv, eps, vscale, vlambda, GHAL, DHAL, SCEXP,
                        SCALPHA, e, de);
         if (rik2 > cut2) {
            real taper, dtaper;
            switch_taper5<do_g>(rik, cut, off, taper, dtaper);
            if CONSTEXPR (do_g)
               de = e * dtaper + de * taper;
            if CONSTEXPR (do_e)
               e = e * taper;
         }


         if CONSTEXPR (do_a)
            if (vscale == -1)
               atomic_add(-1, nev, offset);
         if CONSTEXPR (do_e)
            atomic_add(e, ev, offset);
         if CONSTEXPR (do_g) {
            de *= REAL_RECIP(rik);
            real dedx = de * xr;
            real dedy = de * yr;
            real dedz = de * zr;
            atomic_add(dedx, gxred, i);
            atomic_add(dedy, gyred, i);
            atomic_add(dedz, gzred, i);
            atomic_add(-dedx, gxred, k);
            atomic_add(-dedy, gyred, k);
            atomic_add(-dedz, gzred, k);
            if CONSTEXPR (do_v) {
               real vxx = xr * dedx;
               real vyx = yr * dedx;
               real vzx = zr * dedx;
               real vyy = yr * dedy;
               real vzy = zr * dedy;
               real vzz = zr * dedz;
               atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_ev, offset);
            }
         }
      } // end if (include)
   }
}


template <int USE, evdw_t VDWTYP>
void evdw_cu()
{
   constexpr int do_e = USE & calc::energy;
   constexpr int do_a = USE & calc::analyz;
   constexpr int do_g = USE & calc::grad;
   constexpr int do_v = USE & calc::virial;
   static_assert(do_v ? do_g : true, "");
   static_assert(do_a ? do_e : true, "");

   const auto& st = *vspatial_unit;
   const real cut = switch_cut(switch_vdw);
   const real off = st.cutoff;
   const auto* sp = vspatial_unit.deviceptr();

   auto bufsize = buffer_size();

   if CONSTEXPR (do_g) {
      zero_gradient_async(n, gxred, gyred, gzred);
   }

   if CONSTEXPR (VDWTYP == evdw_t::hal) {
      if (st.niak > 0)
         launch_k1s(nonblk, WARP_SIZE * st.niak, evdw_hal_cu1<USE>, bufsize,
                    nev, ev, vir_ev, gxred, gyred, gzred, TINKER_IMAGE_ARGS,
                    njvdw, jvdw, radmin, epsilon, vlam, vcouple, cut, off, n,
                    st.sorted, st.niak, st.iak, st.lst);
      if (nvexclude_ > 0)
         launch_k1s(nonblk, nvexclude_, evdw_hal_cu2<USE>, bufsize, nev, ev,
                    vir_ev, gxred, gyred, gzred, TINKER_IMAGE_ARGS, njvdw, jvdw,
                    radmin, epsilon, vlam, vcouple, cut, off, xred, yred, zred,
                    nvexclude_, vexclude_, vexclude_scale_);
   }

   if CONSTEXPR (do_g) {
      evdw_resolve_gradient();
   }
}


void evdw_hal_cu(int vers)
{
   if (vers == calc::v0)
      evdw_cu<calc::v0, evdw_t::hal>();
   else if (vers == calc::v1)
      evdw_cu<calc::v1, evdw_t::hal>();
   else if (vers == calc::v3)
      evdw_cu<calc::v3, evdw_t::hal>();
   else if (vers == calc::v4)
      evdw_cu<calc::v4, evdw_t::hal>();
   else if (vers == calc::v5)
      evdw_cu<calc::v5, evdw_t::hal>();
   else if (vers == calc::v6)
      evdw_cu<calc::v6, evdw_t::hal>();
}
TINKER_NAMESPACE_END
