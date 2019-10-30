#pragma once
#include "macro.h"
#include <fftw3.h>


TINKER_NAMESPACE_BEGIN
/// \brief FFT plan.
struct FFTPlan
{
#if TINKER_SINGLE_PRECISION
   using type = fftwf_plan;
#elif TINKER_DOUBLE_PRECISION
   using type = fftw_plan;
#else
   static_assert(false, "");
#endif


   type planf; ///< FFT front plan.
   type planb; ///< FFT back plan.
};
TINKER_NAMESPACE_END
