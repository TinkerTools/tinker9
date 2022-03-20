#pragma once
#include "macro.h"
#include <type_traits>

namespace tinker {
/**
 * \ingroup fft
 * \page fft
 *
 * \warning
 *
 * %PME grid sizes `nfft1`, `nfft2`, and `nfft3` in Tinker are associated with
 * `x`, `y`, and `z` directions, respectively. These three numbers also
 * correspond to a Fortran array `qgrid(2,nfft1,nfft2,nfft3)`, where `nfft1` is
 * the fastest changing dimension and `nfft3` is the slowest changing dimension
 * among the three.
 *
 * <br>
 *
 * The multi-dimensional `FFTW3` Fortran API in Tinker is called as follows,
 * \code{.f}
 *       call dfftw_plan_dft_3d (planf,nfft1,nfft2,nfft3,qgrid,qgrid,ifront,iguess)
 * \endcode
 * which is different from the C API.
 *
 * <br>
 *
 * The difference is described on the [FFTW website]
 * (http://fftw.org/doc/Reversing-array-dimensions.html#Reversing-array-dimensions):
 *
 * > A minor annoyance in calling FFTW from Fortran is that FFTW’s array
 * > dimensions are defined in the C convention (row-major order), while
 * > Fortran’s array dimensions are the opposite convention (column-major
 * > order). See Multi-dimensional Array Format. This is just a bookkeeping
 * > difference, with no effect on performance. The only consequence of this is
 * > that, whenever you create an FFTW plan for a multi-dimensional transform,
 * > you must always reverse the ordering of the dimensions.
 *
 *
 * [`cuFFT` C API]
 * (https://docs.nvidia.com/cuda/cufft/index.html#function-cufftplan3d)
 * is similar to FFTW C API:
 *
 * \code{.cpp}
 * cufftResult
 *    cufftPlan3d(cufftHandle *plan, int nx, int ny, int nz, cufftType type);
 * \endcode
 *
 * > nx: The transform size in the x dimension. This is slowest changing
 * > dimension of a transform (strided in memory).
 * > <br>
 * > nz: The transform size in the z dimension. This is fastest changing
 * > dimension of a transform (contiguous in memory).
 */

struct FFTPlan
{
   template <class T>
   T& self()
   {
      static_assert(std::is_base_of<FFTPlan, T>::value, "");
      return *reinterpret_cast<T*>(this);
   }

   virtual ~FFTPlan() {}
};
}
