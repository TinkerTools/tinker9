#include "acc_add.h"
#include "box.h"
#include "elec.h"
#include "mathfunc.h"
#include "md.h"
#include "pme.h"

TINKER_NAMESPACE_BEGIN
template <int DO_V>
void pme_conv_tmpl(PMEUnit pme_u, virial_buffer gpu_vir)
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

   auto bufsize = buffer_size();
   #pragma acc parallel loop independent async\
               deviceptr(gpu_vir,box,qgrid,bsmod1,bsmod2,bsmod3)
   for (int i = 0; i < ntot; ++i) {
      real recip[3][3];
      recip[0][0] = box->recip[0][0];
      recip[0][1] = box->recip[0][1];
      recip[0][2] = box->recip[0][2];
      recip[1][0] = box->recip[1][0];
      recip[1][1] = box->recip[1][1];
      recip[1][2] = box->recip[1][2];
      recip[2][0] = box->recip[2][0];
      recip[2][1] = box->recip[2][1];
      recip[2][2] = box->recip[2][2];

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

      real h1 = recip[0][0] * r1 + recip[1][0] * r2 + recip[2][0] * r3;
      real h2 = recip[0][1] * r1 + recip[1][1] * r2 + recip[2][1] * r3;
      real h3 = recip[0][2] * r1 + recip[1][2] * r2 + recip[2][2] * r3;
      real hsq = h1 * h1 + h2 * h2 + h3 * h3;

      real gridx = qgrid[i][0];
      real gridy = qgrid[i][1];
      real term = -pterm * hsq;
      real expterm = 0;
      if (term > -50) {
         // TODO: if .not. use_bounds; if octahedron; 2/hsq
         real denom =
            hsq * pi * box->volbox * bsmod1[k1] * bsmod2[k2] * bsmod3[k3];
         expterm = REAL_EXP(term) / denom;

         if CONSTEXPR (DO_V) {
            real struc2 = gridx * gridx + gridy * gridy;
            real eterm = 0.5f * f * expterm * struc2;
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
      }

      // complete the transformation of the PME grid
      qgrid[i][0] = gridx * expterm;
      qgrid[i][1] = gridy * expterm;
   }
}

void pme_conv0(PMEUnit pme_u)
{
   pme_conv_tmpl<0>(pme_u, nullptr);
}

void pme_conv1(PMEUnit pme_u, virial_buffer gpu_vir)
{
   pme_conv_tmpl<1>(pme_u, gpu_vir);
}
TINKER_NAMESPACE_END
