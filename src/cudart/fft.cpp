#include "fft.h"
#include "cudalib.h"
#include "error.h"
#include "pme.h"
#include <cufft.h>


TINKER_NAMESPACE_BEGIN
struct FFTPlanCUFFT : public FFTPlan
{
   cufftHandle h;
};


void fft_data(rc_op op)
{
   if (op & rc_dealloc) {
      int idx = 0;
      while (idx < FFTPlanUnit::size()) {
         FFTPlanUnit u = idx;
         cufftDestroy(u->self<FFTPlanCUFFT>().h);
         ++idx;
      }
      FFTPlanUnit::clear();
   }

   if (op & rc_alloc) {
      assert(FFTPlanUnit::size() == 0);

      const size_t size = PMEUnit::size();
      FFTPlanUnit::resize<FFTPlanCUFFT>(size);
   }

   if (op & rc_init) {
#if TINKER_SINGLE_PRECISION
      const cufftType typ = CUFFT_C2C;
#elif TINKER_DOUBLE_PRECISION
      const cufftType typ = CUFFT_Z2Z;
#else
      static_assert(false, "");
#endif

      int idx = 0;
      while (idx < FFTPlanUnit::size()) {
         FFTPlanUnit plan_u = idx;
         PMEUnit pme_u = idx;
         auto& iplan = plan_u->self<FFTPlanCUFFT>().h;
         auto& st = *pme_u;

         check_rt(cufftPlan3d(&iplan, st.nfft1, st.nfft2, st.nfft3, typ));
         check_rt(cufftSetStream(iplan, nonblk));
         ++idx;
      }
   }
}


void fftfront(PMEUnit pme_u)
{
   FFTPlanUnit iplan_u = static_cast<int>(pme_u);
   auto& iplan = iplan_u->self<FFTPlanCUFFT>().h;
   auto& st = *pme_u;

#if TINKER_SINGLE_PRECISION
   cufftExecC2C(iplan, reinterpret_cast<cufftComplex*>(st.qgrid),
                reinterpret_cast<cufftComplex*>(st.qgrid), CUFFT_FORWARD);
#elif TINKER_DOUBLE_PRECISION
   cufftExecZ2Z(iplan, reinterpret_cast<cufftDoubleComplex*>(st.qgrid),
                reinterpret_cast<cufftDoubleComplex*>(st.qgrid), CUFFT_FORWARD);
#else
   static_assert(false, "");
#endif
}


void fftback(PMEUnit pme_u)
{
   FFTPlanUnit iplan_u = static_cast<int>(pme_u);
   auto& iplan = iplan_u->self<FFTPlanCUFFT>().h;
   auto& st = *pme_u;

#if TINKER_SINGLE_PRECISION
   cufftExecC2C(iplan, reinterpret_cast<cufftComplex*>(st.qgrid),
                reinterpret_cast<cufftComplex*>(st.qgrid), CUFFT_INVERSE);
#elif TINKER_DOUBLE_PRECISION
   cufftExecZ2Z(iplan, reinterpret_cast<cufftDoubleComplex*>(st.qgrid),
                reinterpret_cast<cufftDoubleComplex*>(st.qgrid), CUFFT_INVERSE);
#else
   static_assert(false, "");
#endif
}
TINKER_NAMESPACE_END
