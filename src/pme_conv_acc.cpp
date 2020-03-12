#include "add.h"
#include "box.h"
#include "elec.h"
#include "mathfunc.h"
#include "md.h"
#include "pmestuf.h"


TINKER_NAMESPACE_BEGIN
template <bool DO_E, bool DO_V>
void pme_conv_acc1(PMEUnit pme_u, energy_buffer gpu_e, virial_buffer gpu_vir)
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


   const real f = electric / dielec;
   real aewald = st.aewald;
   real pterm = pi / aewald;
   pterm *= pterm;
   real box_volume = volbox();


   auto bufsize = buffer_size();
   #pragma acc parallel loop independent async\
               deviceptr(gpu_e,gpu_vir,qgrid,bsmod1,bsmod2,bsmod3)
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
      real term = -pterm * hsq;
      real expterm = 0;
      if (term > -50) {
         // TODO: if .not. use_bounds; if octahedron; 2/hsq
         real denom =
            hsq * pi * box_volume * bsmod1[k1] * bsmod2[k2] * bsmod3[k3];
         expterm = REAL_EXP(term) / denom;


         if CONSTEXPR (DO_E || DO_V) {
            real struc2 = gridx * gridx + gridy * gridy;
            real eterm = 0.5f * f * expterm * struc2;
            if (DO_E) {
               atomic_add(eterm, gpu_e, i & (bufsize - 1));
            }
            if (DO_V) {
               real vterm = (2 / hsq) * (1 - term) * eterm;

               real vxx = (h1 * h1 * vterm - eterm);
               real vxy = h1 * h2 * vterm;
               real vxz = h1 * h3 * vterm;
               real vyy = (h2 * h2 * vterm - eterm);
               real vyz = h2 * h3 * vterm;
               real vzz = (h3 * h3 * vterm - eterm);
               atomic_add(vxx, vxy, vxz, vyy, vyz, vzz, gpu_vir,
                          i & (bufsize - 1));
            }
         } // end if (e or v)
      }


      // complete the transformation of the PME grid
      qgrid[i][0] = gridx * expterm;
      qgrid[i][1] = gridy * expterm;
   }
}


void pme_conv_acc(PMEUnit pme_u, energy_buffer gpu_e, virial_buffer gpu_vir)
{
   if (gpu_vir == nullptr) {
      if (gpu_e == nullptr) {
         pme_conv_acc1<false, false>(pme_u, nullptr, nullptr);
      } else {
         pme_conv_acc1<true, false>(pme_u, gpu_e, nullptr);
      }
   } else {
      if (gpu_e == nullptr) {
         pme_conv_acc1<false, true>(pme_u, nullptr, gpu_vir);
      } else {
         pme_conv_acc1<true, true>(pme_u, gpu_e, gpu_vir);
      }
   }
}
TINKER_NAMESPACE_END
