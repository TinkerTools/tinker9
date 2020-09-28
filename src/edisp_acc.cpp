#include "add.h"
#include "box.h"
#include "edisp.h"
#include "md.h"
#include "pmestuf.h"


namespace tinker {
template <bool DO_E, bool DO_V>
void disp_pme_conv_acc1(PMEUnit pme_u, energy_buffer gpu_e, virial_buffer gpu_v)
{
   auto& st = *pme_u;
   real(*restrict qgrid)[2] = reinterpret_cast<real(*)[2]>(st.qgrid);
   const real* bsmod1 = st.bsmod1;
   const real* bsmod2 = st.bsmod2;
   const real* bsmod3 = st.bsmod3;


   const int nfft1 = st.nfft1;
   const int nfft2 = st.nfft2;
   const int nfft3 = st.nfft3;
   const int nff = nfft1 * nfft2;
   const int ntot = nfft1 * nfft2 * nfft3;


   const real aewald = st.aewald;
   const real bfac = M_PI / aewald;
   const real fac1 = 2 * std::pow(M_PI, 3.5);
   const real fac2 = aewald * aewald * aewald;
   const real fac3 = -2 * aewald * M_PI * M_PI;
   const real vbox = volbox();
   const real denom0 = 6 * vbox / std::pow(M_PI, 1.5);


   size_t bufsize = buffer_size();
   #pragma acc parallel loop independent async\
               deviceptr(gpu_e,gpu_v,qgrid,bsmod1,bsmod2,bsmod3)
   for (int i = 0; i < ntot; ++i) {
      if (i == 0) {
         qgrid[0][0] = 0;
         qgrid[0][1] = 0;
         continue;
      }


      int k3 = i / nff;
      int j = i - k3 * nff;
      int k2 = j / nfft1;
      int k1 = j - k2 * nfft1;


      int r1 = (k1 < (nfft1 + 1) / 2) ? k1 : (k1 - nfft1);
      int r2 = (k2 < (nfft2 + 1) / 2) ? k2 : (k2 - nfft2);
      int r3 = (k3 < (nfft3 + 1) / 2) ? k3 : (k3 - nfft3);


      real h1 = recipa.x * r1 + recipb.x * r2 + recipc.x * r3;
      real h2 = recipa.y * r1 + recipb.y * r2 + recipc.y * r3;
      real h3 = recipa.z * r1 + recipb.z * r2 + recipc.z * r3;
      real hsq = h1 * h1 + h2 * h2 + h3 * h3;


      real gridx = qgrid[i][0];
      real gridy = qgrid[i][1];
      real h = REAL_SQRT(hsq);
      real b = h * bfac;
      real hhh = h * hsq;
      real term = -hsq * bfac * bfac;
      real eterm = 0;
      real denom = denom0 * bsmod1[k1] * bsmod2[k2] * bsmod3[k3];
      if (term > -50) {
         real expterm = REAL_EXP(term);
         real erfcterm = REAL_ERFC(b);
         if (box_shape == UNBOUND_BOX) {
            real coef = (1 - REAL_COS(pi * lvec1.x * REAL_SQRT(hsq)));
            expterm *= coef;
            erfcterm *= coef;
         } else if (box_shape == OCT_BOX) {
            if ((k1 + k2 + k3) & 1) {
               expterm = 0;
               erfcterm = 0;
            } // end if ((k1 + k2 + k3) % 2 != 0)
         }


         real struc2 = gridx * gridx + gridy * gridy;
         eterm = (-fac1 * erfcterm * hhh - expterm * (fac2 + fac3 * hsq)) *
            REAL_RECIP(denom);
         real e = eterm * struc2;
         if CONSTEXPR (DO_E) {
            atomic_add(e, gpu_e, i & (bufsize - 1));
         }
         if CONSTEXPR (DO_V) {
            real vterm = 3 * (fac1 * erfcterm * h + fac3 * expterm) * struc2 *
               REAL_RECIP(denom);
            real vxx = (h1 * h1 * vterm - e);
            real vxy = h1 * h2 * vterm;
            real vxz = h1 * h3 * vterm;
            real vyy = (h2 * h2 * vterm - e);
            real vyz = h2 * h3 * vterm;
            real vzz = (h3 * h3 * vterm - e);
            atomic_add(vxx, vxy, vxz, vyy, vyz, vzz, gpu_v, i & (bufsize - 1));
         }
      }


      qgrid[i][0] = eterm * gridx;
      qgrid[i][1] = eterm * gridy;
   }
}


void disp_pme_conv_acc(int vers)
{
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   PMEUnit u = dpme_unit;


   if (do_e && do_v)
      disp_pme_conv_acc1<true, true>(u, edsp, vir_edsp);
   else if (do_e && !do_v)
      disp_pme_conv_acc1<true, false>(u, edsp, nullptr);
   else if (!do_e && do_v)
      disp_pme_conv_acc1<false, true>(u, nullptr, vir_edsp);
   else if (!do_e && !do_v)
      disp_pme_conv_acc1<false, false>(u, nullptr, nullptr);
}
}
