#ifndef TINKER_FFT_FFTW_H_
#define TINKER_FFT_FFTW_H_

#include "macro.h"
#include <fftw3.h>

TINKER_NAMESPACE_BEGIN
/// @brief
/// FFT plan
struct FFTPlan {
#if TINKER_SINGLE_PRECISION
  typedef fftwf_plan type;
#elif TINKER_DOUBLE_PRECISION
  typedef fftw_plan type;
#else
  static_assert(false, "");
#endif

  type planf; ///< FFT front plan
  type planb; ///< FFT back plan
};
TINKER_NAMESPACE_END

#endif
