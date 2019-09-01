#ifndef TINKER_FFT_H_
#define TINKER_FFT_H_

#ifdef TINKER_CUDART
#  include "fft_cufft.h"
#endif

#ifdef TINKER_HOST
#  include "fft_fftw.h"
#endif

#endif
