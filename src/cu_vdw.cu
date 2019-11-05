#include "cu_add.h"
#include "cu_launch.h"
#include "e_vdw.h"
#include "md.h"
#include "seq_image.h"
#include "seq_pair_hal.h"
#include "seq_switch.h"
#include "spatial.h"


/**
 * Overheads:
 *    - Different vcouple methods.
 *    - PBC type in image().
 *    - (If not hard-coded) ghal, dhal, scexp, scalpha.
 */


TINKER_NAMESPACE_BEGIN
#define HAL_ARGS                                                               \
   size_t bufsize, count_buffer restrict nev, energy_buffer restrict ev,       \
      virial_buffer vir_ev, real *restrict gxred, real *restrict gyred,        \
      real *restrict gzred, const Box *restrict box, int njvdw,                \
      const int *restrict jvdw, const real *restrict radmin,                   \
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
__global__
void evdw_hal_cu1(HAL_ARGS, const Spatial* restrict sp)
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
   MAYBE_UNUSED real gxi, gyi, gzi;
   MAYBE_UNUSED __shared__ real gxk[BLOCK_DIM], gyk[BLOCK_DIM], gzk[BLOCK_DIM];
   MAYBE_UNUSED real vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz;


   const real cut2 = cut * cut;
   const real off2 = off * off;
   const int n = sp->n;
   const int niak = sp->niak;
   const auto* restrict sorted = sp->sorted;
   const auto* restrict iak = sp->iak;
   const auto* restrict lst = sp->lst;
   for (int iw = iwarp; iw < niak; iw += nwarp) {
      /**
       * ATOMI may exceed n-1 for the last AtomBlock: set it to n-1 if that is
       * the case. On the other hand, ATOMK has been stored in the list and is
       * guaranteed that it is not greater than n-1, therefore,
       * `IF (ATOMI < ATOMK)` will exclude all of the "overflowed" atom pairs.
       */
      int atomi; // sorted atom number i
      atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
      real xi = sorted[atomi].x;
      real yi = sorted[atomi].y;
      real zi = sorted[atomi].z;
      int i = sorted[atomi].unsorted;
      int it = jvdw[i];
      real lam1 = vlam[i];


      int shatomk; // sorted atom number k stored in ilane
      shatomk = lst[iw * WARP_SIZE + ilane];
      real shx = sorted[shatomk].x;
      real shy = sorted[shatomk].y;
      real shz = sorted[shatomk].z;
      int shk = sorted[shatomk].unsorted;
      int shkt = jvdw[shk];
      real shlam = vlam[shk];


      if_constexpr(do_a) ctl = 0;
      if_constexpr(do_e) etl = 0;
      if_constexpr(do_g)
      {
         gxi = 0;
         gyi = 0;
         gzi = 0;
         gxk[threadIdx.x] = 0;
         gyk[threadIdx.x] = 0;
         gzk[threadIdx.x] = 0;
      }
      if_constexpr(do_v)
      {
         vtlxx = 0;
         vtlyx = 0;
         vtlzx = 0;
         vtlyy = 0;
         vtlzy = 0;
         vtlzz = 0;
      }


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
         real xr = xi - __shfl_sync(ALL_LANES, shx, srclane);
         real yr = yi - __shfl_sync(ALL_LANES, shy, srclane);
         real zr = zi - __shfl_sync(ALL_LANES, shz, srclane);
         int k = __shfl_sync(ALL_LANES, shk, srclane);
         int kt = __shfl_sync(ALL_LANES, shkt, srclane);
         real vlambda = __shfl_sync(ALL_LANES, shlam, srclane);
         if (vcouple == evdw_t::decouple) {
            vlambda = (lam1 == vlambda ? 1 : (lam1 < vlambda ? lam1 : vlambda));
         } else if (vcouple == evdw_t::annihilate) {
            vlambda = (lam1 < vlambda ? lam1 : vlambda);
         }


         MAYBE_UNUSED real dedx = 0, dedy = 0, dedz = 0;
         image(xr, yr, zr, box);
         real rik2 = xr * xr + yr * yr + zr * zr;
         if (atomi < atomk && rik2 <= off2) {
            real rik = REAL_SQRT(rik2);
            real rv = radmin[it * njvdw + kt];
            real eps = epsilon[it * njvdw + kt];


            MAYBE_UNUSED real e, de;
            pair_hal<do_g>(rik, rv, eps, 1, vlambda, GHAL, DHAL, SCEXP, SCALPHA,
                           e, de);
            if (rik2 > cut2) {
               real taper, dtaper;
               switch_taper5<do_g>(rik, cut, off, taper, dtaper);
               if_constexpr(do_g) de = e * dtaper + de * taper;
               if_constexpr(do_e) e = e * taper;
            }


            if_constexpr(do_a) ctl += 1;
            if_constexpr(do_e) etl += e;
            if_constexpr(do_g)
            {
               de *= REAL_RECIP(rik);
               dedx = de * xr;
               dedy = de * yr;
               dedz = de * zr;


               if_constexpr(do_v)
               {
                  vtlxx += xr * dedx;
                  vtlyx += yr * dedx;
                  vtlzx += zr * dedx;
                  vtlyy += yr * dedy;
                  vtlzy += zr * dedy;
                  vtlzz += zr * dedz;
               } // end if (do_v)
            }    // end if (do_g)
         }


         if_constexpr(do_g)
         {
            gxi += dedx;
            gyi += dedy;
            gzi += dedz;
            gxk[srclane + (threadIdx.x - ilane)] -= dedx;
            gyk[srclane + (threadIdx.x - ilane)] -= dedy;
            gzk[srclane + (threadIdx.x - ilane)] -= dedz;
         }
      }


      if_constexpr(do_a) atomic_add_value(ctl, nev, offset);
      if_constexpr(do_e) atomic_add_value(etl, ev, offset);
      if_constexpr(do_g)
      {
         atomic_add_value(gxi, gxred, i);
         atomic_add_value(gyi, gyred, i);
         atomic_add_value(gzi, gzred, i);
         atomic_add_value(gxk[threadIdx.x], gxred, shk);
         atomic_add_value(gyk[threadIdx.x], gyred, shk);
         atomic_add_value(gzk[threadIdx.x], gzred, shk);
      }
      if_constexpr(do_v) atomic_add_value(vtlxx, vtlyx, vtlzx, vtlyy, vtlzy,
                                          vtlzz, vir_ev, offset);
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
      if (vcouple == evdw_t::decouple) {
         vlambda = (lam1 == vlambda ? 1 : (lam1 < vlambda ? lam1 : vlambda));
      } else if (vcouple == evdw_t::annihilate) {
         vlambda = (lam1 < vlambda ? lam1 : vlambda);
      }

      image(xr, yr, zr, box);
      real rik2 = xr * xr + yr * yr + zr * zr;
      if (rik2 <= off2) {
         real rik = REAL_SQRT(rik2);
         real rv = radmin[it * njvdw + kt];
         real eps = epsilon[it * njvdw + kt];

         MAYBE_UNUSED real e, de;
         pair_hal<do_g>(rik, rv, eps, vscale, vlambda, GHAL, DHAL, SCEXP,
                        SCALPHA, e, de);

         if (rik2 > cut2) {
            real taper, dtaper;
            switch_taper5<do_g>(rik, cut, off, taper, dtaper);
            if_constexpr(do_g) de = e * dtaper + de * taper;
            if_constexpr(do_e) e = e * taper;
         }

         if_constexpr(do_a) if (vscale == -1) atomic_add_value(-1, nev, offset);
         if_constexpr(do_e) atomic_add_value(e, ev, offset);
         if_constexpr(do_g)
         {
            de *= REAL_RECIP(rik);
            real dedx = de * xr;
            real dedy = de * yr;
            real dedz = de * zr;

            atomic_add_value(dedx, gxred, i);
            atomic_add_value(dedy, gyred, i);
            atomic_add_value(dedz, gzred, i);
            atomic_add_value(-dedx, gxred, k);
            atomic_add_value(-dedy, gyred, k);
            atomic_add_value(-dedz, gzred, k);

            if_constexpr(do_v)
            {
               real vxx = xr * dedx;
               real vyx = yr * dedx;
               real vzx = zr * dedx;
               real vyy = yr * dedy;
               real vzy = zr * dedy;
               real vzz = zr * dedz;
               atomic_add_value(vxx, vyx, vzx, vyy, vzy, vzz, vir_ev, offset);
            } // end if (do_v)
         }    // end if (do_g)
      }
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

   const real cut = switch_cut(switch_vdw);
   const real off = switch_off(switch_vdw);
   const auto& st = *vspatial_unit;
   const auto* sp = vspatial_unit.deviceptr();

   auto bufsize = buffer_size();

   if_constexpr(do_g)
   {
      device_array::zero(n, gxred, gyred, gzred);
   }

   if_constexpr(VDWTYP == evdw_t::hal)
   {
      launch_kernel1(WARP_SIZE * st.niak, evdw_hal_cu1<USE>, bufsize, nev, ev,
                     vir_ev, gxred, gyred, gzred, box, njvdw, jvdw, radmin,
                     epsilon, vlam, vcouple, cut, off, sp);
      if (nvexclude_ > 0)
         launch_kernel1(nvexclude_, evdw_hal_cu2<USE>, bufsize, nev, ev, vir_ev,
                        gxred, gyred, gzred, box, njvdw, jvdw, radmin, epsilon,
                        vlam, vcouple, cut, off, xred, yred, zred, nvexclude_,
                        vexclude_, vexclude_scale_);
   }

   if_constexpr(do_g)
   {
      evdw_resolve_gradient();
   }
}


void evdw_hal_cu(int vers)
{
   evdw_reduce_xyz();
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
