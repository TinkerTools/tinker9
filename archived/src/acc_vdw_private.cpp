#include "spatial.h"
#ifdef _OPENACC
#include <openacc.h>
#endif


#define GHAL (real)0.12
#define DHAL (real)0.07
#define SCEXP 5
#define SCALPHA (real)0.7
template <int USE, evdw_t VDWTYP>
void evdw_tmpl()
{
   constexpr int do_e = USE & calc::energy;
   constexpr int do_a = USE & calc::analyz;
   constexpr int do_g = USE & calc::grad;
   constexpr int do_v = USE & calc::virial;
   if CONSTEXPR (do_g)
      zero_gradient(false, n, gxred, gyred, gzred);

   const auto& st = *vspatial_unit;
   const auto* sorted = st.sorted;
   const auto* iak = st.iak;
   const auto* lst = st.lst;
   int niak = st.niak;

   const real cut = switch_cut(switch_vdw);
   const real off = st.cutoff;
   const real cut2 = cut * cut;
   const real off2 = off * off;
   auto bufsize = buffer_size();

   constexpr int BLKS = WARP_SIZE;
   real gxk[BLKS], gyk[BLKS], gzk[BLKS];
   real shx[BLKS], shy[BLKS], shz[BLKS];
   real shlam[BLKS];
   int shatomk[BLKS], shk[BLKS], shkt[BLKS];

   //*
   #pragma acc parallel loop async gang\
               deviceptr(sorted,iak,lst,xred,yred,zred,gxred,gyred,gzred,\
               jvdw,radmin,epsilon,vlam,nev,ev,vir_ev) \
               private(shx[BLKS],shy[BLKS],shz[BLKS],\
               shlam[BLKS],shatomk[BLKS],shk[BLKS],shkt[BLKS],gxk[BLKS],gyk[BLKS],gzk[BLKS])\
               vector_length(WARP_SIZE)
   for (int iw = 0; iw < niak; ++iw) {
      // load data into shared memory
      #pragma acc loop vector
      for (int t = 0; t < BLKS; ++t) {
         if CONSTEXPR (do_g) {
            gxk[t] = 0;
            gyk[t] = 0;
            gzk[t] = 0;
         }
         int ilane = t;
         int atomk = lst[iw * WARP_SIZE + ilane];
         shx[t] = sorted[atomk].x;
         shy[t] = sorted[atomk].y;
         shz[t] = sorted[atomk].z;
         int k = sorted[atomk].unsorted;
         shatomk[t] = atomk;
         shk[t] = k;
         shkt[t] = jvdw[k];
         shlam[t] = vlam[k];
      }
      // implicit __syncthreads();

      #pragma acc loop vector
      for (int t = 0; t < BLKS; ++t) {
         MAYBE_UNUSED int ctl;
         MAYBE_UNUSED real etl;
         MAYBE_UNUSED real gxi, gyi, gzi;
         MAYBE_UNUSED real vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz;
         if CONSTEXPR (do_e)
            ctl = 0;
         if CONSTEXPR (do_e)
            etl = 0;
         if CONSTEXPR (do_g) {
            gxi = 0;
            gyi = 0;
            gzi = 0;
         }
         if CONSTEXPR (do_v) {
            vtlxx = 0;
            vtlyx = 0;
            vtlzx = 0;
            vtlyy = 0;
            vtlzy = 0;
            vtlzz = 0;
         }

         int offset = (t + __pgi_gangidx() * BLKS) & (bufsize - 1);
         int ilane = t;
         int ai0 = iak[iw] * WARP_SIZE + ilane;
         int atomi = ai0 < n ? ai0 : (n - 1);
         real xi = sorted[atomi].x;
         real yi = sorted[atomi].y;
         real zi = sorted[atomi].z;
         int i = sorted[atomi].unsorted;
         int it = jvdw[i];
         real lam1 = vlam[i];

         #pragma acc loop seq
         for (int j = 0; j < WARP_SIZE; ++j) {
            int srclane = (ilane + j) & 31;
            int klane = srclane;
            int atomk = shatomk[klane];
            real xr = xi - shx[klane];
            real yr = yi - shy[klane];
            real zr = zi - shz[klane];
            int kt = shkt[klane];
            real vlambda = shlam[klane];

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
               pair_hal<do_g>(rik, rv, eps, 1, vlambda, GHAL, DHAL, SCEXP,
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
                  ctl += 1;
               if CONSTEXPR (do_e)
                  etl += e;
               if CONSTEXPR (do_g) {
                  de *= REAL_RECIP(rik);
                  real dedx = de * xr;
                  real dedy = de * yr;
                  real dedz = de * zr;
                  gxi += dedx;
                  gyi += dedy;
                  gzi += dedz;
                  gxk[klane] -= dedx;
                  gyk[klane] -= dedy;
                  gzk[klane] -= dedz;
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
         }

         if CONSTEXPR (do_a)
            atomic_add(ctl, nev, offset);
         if CONSTEXPR (do_e)
            atomic_add(etl, ev, offset);
         if CONSTEXPR (do_g) {
            atomic_add(gxi, gxred, i);
            atomic_add(gyi, gyred, i);
            atomic_add(gzi, gzred, i);
            atomic_add(gxk[t], gxred, shk[t]);
            atomic_add(gyk[t], gyred, shk[t]);
            atomic_add(gzk[t], gzred, shk[t]);
         }
         if CONSTEXPR (do_v)
            atomic_add(vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz, vir_ev,
                       offset);
      }
   }
   //*/

   //*
   #pragma acc parallel async deviceptr(xred,yred,zred,gxred,gyred,gzred,\
               jvdw,radmin,epsilon,vlam,\
               nev,ev,vir_ev,vexclude_,vexclude_scale_) 
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

      real rik2 = image2(xr, yr, zr);
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
   // */

   if CONSTEXPR (do_g)
      evdw_resolve_gradient();
}
