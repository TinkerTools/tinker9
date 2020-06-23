#include "add.h"
#include "box.h"
#include "edisp.h"
#include "pmestuf.h"


namespace tinker {
template <bool DO_E, bool DO_V>
void disppme_conv_acc1(PMEUnit pme_u, energy_buffer gpu_e, virial_buffer gpu_v)
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


      real h = REAL_SQRT(hsq);
      real b = h * bfac;
      real hhh = h * hsq;
      real term = -hsq * bfac * bfac;
      real expterm = 0;
      real erfcterm = REAL_ERFC(b);
      real denom = denom0 * bsmod1[k1] * bsmod2[k2] * bsmod3[k3];
      if (term > -50) {
      }
   }
}


template <class Ver>
void edisp_ewald_recip_self_acc1()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;


   PMEUnit u = dpme_unit;
   grid_disp(u, csix);
   fftfront(u);


   if CONSTEXPR (do_e) {
      const real aewald = u->aewald;
      real term = aewald * aewald * aewald;
      term *= term;
      term /= 12;
      size_t bufsize = buffer_size();
      #pragma acc parallel loop independent async\
                  deviceptr(csix,edsp,ndisp)
      for (int i = 0; i < n; ++i) {
         int offset = i & (bufsize - 1);
         real icsix = csix[i];
         atomic_add(term * icsix * icsix, edsp, offset);
         if CONSTEXPR (do_a)
            atomic_add(1, ndisp, offset);
      }
   }
}


void edisp_ewald_recip_self_acc(int vers) {}
}
