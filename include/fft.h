#ifndef TINKER_FFT_H_
#define TINKER_FFT_H_

#include "macro.h"

#if TINKER_CUDART
#  include "fft_cufft.h"
#endif

#if TINKER_HOST
#  include "fft_fftw.h"
#endif

#endif
