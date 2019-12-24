#include "add.h"
#include "e_vdw.h"
#include "gpu_card.h"
#include "md.h"
#include "nblist.h"
#include "seq_image.h"
#include "seq_pair_hal.h"
#include "seq_switch.h"

TINKER_NAMESPACE_BEGIN
void evdw_reduce_xyz()
{
   #pragma acc parallel loop independent\
               deviceptr(x,y,z,ired,kred,xred,yred,zred)
   for (int i = 0; i < n; ++i) {
      int iv = ired[i];
      real rdn = kred[i];
      xred[i] = rdn * (x[i] - x[iv]) + x[iv];
      yred[i] = rdn * (y[i] - y[iv]) + y[iv];
      zred[i] = rdn * (z[i] - z[iv]) + z[iv];
   }
}

void evdw_resolve_gradient()
{
   #pragma acc parallel loop independent async\
               deviceptr(ired,kred,gxred,gyred,gzred,gx,gy,gz)
   for (int ii = 0; ii < n; ++ii) {
      int iv = ired[ii];
      real fx = gxred[ii];
      real fy = gyred[ii];
      real fz = gzred[ii];
      if (ii == iv) {
         atomic_add(fx, gx, ii);
         atomic_add(fy, gy, ii);
         atomic_add(fz, gz, ii);
      } else {
         real redii = kred[ii];
         real rediv = 1 - redii;
         atomic_add(fx * redii, gx, ii);
         atomic_add(fy * redii, gy, ii);
         atomic_add(fz * redii, gz, ii);
         atomic_add(fx * rediv, gx, iv);
         atomic_add(fy * rediv, gy, iv);
         atomic_add(fz * rediv, gz, iv);
      }
   }
}

template <int USE, evdw_t VDWTYP>
void evdw_tmpl()
{
   constexpr int do_e = USE & calc::energy;
   constexpr int do_a = USE & calc::analyz;
   constexpr int do_g = USE & calc::grad;
   constexpr int do_v = USE & calc::virial;
   static_assert(do_v ? do_g : true, "");
   static_assert(do_a ? do_e : true, "");

   const real cut = switch_cut(switch_vdw);
   const real off = vlist_unit->cutoff;
   const real cut2 = cut * cut;
   const real off2 = off * off;
   const int maxnlst = vlist_unit->maxnlst;
   const auto* vlst = vlist_unit.deviceptr();

   auto bufsize = buffer_size();

   if CONSTEXPR (do_g)
      zero_gradient(n, gxred, gyred, gzred);

#define DEVICE_PTRS_                                                           \
   xred, yred, zred, gxred, gyred, gzred, box, jvdw, radmin, epsilon, vlam,    \
      nev, ev, vir_ev

   MAYBE_UNUSED int GRID_DIM = get_grid_size(BLOCK_DIM);
   #pragma acc parallel num_gangs(GRID_DIM) vector_length(BLOCK_DIM)\
               deviceptr(DEVICE_PTRS_,vlst)
   #pragma acc loop gang independent
   for (int i = 0; i < n; ++i) {
      int it = jvdw[i];
      real xi = xred[i];
      real yi = yred[i];
      real zi = zred[i];
      real lam1 = vlam[i];
      MAYBE_UNUSED real gxi = 0, gyi = 0, gzi = 0;

      int nvlsti = vlst->nlst[i];
      #pragma acc loop vector independent reduction(+:gxi,gyi,gzi)
      for (int kk = 0; kk < nvlsti; ++kk) {
         int offset = (kk + i * n) & (bufsize - 1);
         int k = vlst->lst[i * maxnlst + kk];
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
            if CONSTEXPR (VDWTYP == evdw_t::hal)
               pair_hal<do_g>(rik, rv, eps, 1, vlambda,   //
                              ghal, dhal, scexp, scalpha, //
                              e, de);

            if (rik2 > cut2) {
               real taper, dtaper;
               switch_taper5<do_g>(rik, cut, off, taper, dtaper);
               if CONSTEXPR (do_g)
                  de = e * dtaper + de * taper;
               if CONSTEXPR (do_e)
                  e = e * taper;
            }

            // Increment the energy, gradient, and virial.

            if CONSTEXPR (do_a)
               atomic_add(1, nev, offset);
            if CONSTEXPR (do_e)
               atomic_add(e, ev, offset);
            if CONSTEXPR (do_g) {
               de *= REAL_RECIP(rik);
               real dedx = de * xr;
               real dedy = de * yr;
               real dedz = de * zr;

               gxi += dedx;
               gyi += dedy;
               gzi += dedz;
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
               } // end if (do_v)
            }    // end if (do_g)
         }
      } // end for (int kk)

      if CONSTEXPR (do_g) {
         atomic_add(gxi, gxred, i);
         atomic_add(gyi, gyred, i);
         atomic_add(gzi, gzred, i);
      }
   } // end for (int i)

   #pragma acc parallel deviceptr(DEVICE_PTRS_,vexclude_,vexclude_scale_)
   #pragma acc loop independent
   for (int ii = 0; ii < nvexclude_; ++ii) {
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
         if CONSTEXPR (VDWTYP == evdw_t::hal)
            pair_hal<do_g>(rik, rv, eps, vscale, vlambda, //
                           ghal, dhal, scexp, scalpha,    //
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
            } // end if (do_v)
         }    // end if (do_g)
      }
   } // end for (int ii)

   if CONSTEXPR (do_g)
      evdw_resolve_gradient();
}

#define TINKER_EVDW_ACC(typ)                                                   \
   void evdw_##typ##_acc(int vers)                                             \
   {                                                                           \
      if (vers == calc::v0)                                                    \
         evdw_tmpl<calc::v0, evdw_t::typ>();                                   \
      else if (vers == calc::v1)                                               \
         evdw_tmpl<calc::v1, evdw_t::typ>();                                   \
      else if (vers == calc::v3)                                               \
         evdw_tmpl<calc::v3, evdw_t::typ>();                                   \
      else if (vers == calc::v4)                                               \
         evdw_tmpl<calc::v4, evdw_t::typ>();                                   \
      else if (vers == calc::v5)                                               \
         evdw_tmpl<calc::v5, evdw_t::typ>();                                   \
      else if (vers == calc::v6)                                               \
         evdw_tmpl<calc::v6, evdw_t::typ>();                                   \
   }
TINKER_EVDW_ACC(lj);
TINKER_EVDW_ACC(buck);
TINKER_EVDW_ACC(mm3hb);
TINKER_EVDW_ACC(hal);
TINKER_EVDW_ACC(gauss);
TINKER_NAMESPACE_END
