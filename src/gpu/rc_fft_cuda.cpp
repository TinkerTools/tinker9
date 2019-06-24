
#ifndef TINKER_HOSTONLY

#  include "gpu/decl_pme.h"
#  include "gpu/rc.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
std::vector<fft_plan_t>& fft_plans();
}
TINKER_NAMESPACE_END

//======================================================================

TINKER_NAMESPACE_BEGIN
namespace gpu {
void fft_data(rc_t rc) {
  if (rc & rc_dealloc) {
    int idx = 0;
    while (idx < fft_plans().size()) {
      cufftDestroy(fft_plans()[idx]);
      ++idx;
    }
    fft_plans().clear();
  }

  if (rc & rc_alloc) {
    assert(fft_plans().size() == 0);

    const size_t size = detail_::pme_objs().size();
    fft_plans().resize(size, fft_plan_t());
  }

  if (rc & rc_copyin) {
#  if defined(TINKER_GPU_SINGLE)
    const cufftType typ = CUFFT_C2C;
#  elif defined(TINKER_GPU_DOUBLE)
    const cufftType typ = CUFFT_Z2Z;
#  else
    static_assert(false, "");
#  endif

    int idx = 0;
    for (fft_plan_t& iplan : fft_plans()) {
      auto& st = detail_::pme_objs()[idx];
      check_cudart(cufftPlan3d(&iplan, st.nfft1, st.nfft2, st.nfft3, typ),
                   cufftResult, CUFFT_SUCCESS);
      ++idx;
    }
  }
}

void fftfront(int pme_unit) {
  fft_plan_t iplan = fft_plans()[pme_unit];
  auto& st = detail_::pme_objs()[pme_unit];

#  if defined(TINKER_GPU_SINGLE)
  cufftExecC2C(iplan, reinterpret_cast<cufftComplex*>(st.qgrid),
               reinterpret_cast<cufftComplex*>(st.qgrid), CUFFT_FORWARD);
#  elif defined(TINKER_GPU_DOUBLE)
  cufftExecZ2Z(iplan, reinterpret_cast<cufftDoubleComplex*>(st.qgrid),
               reinterpret_cast<cufftDoubleComplex*>(st.qgrid), CUFFT_FORWARD);
#  else
  static_assert(false, "");
#  endif
}

void fftback(int pme_unit) {
  fft_plan_t iplan = fft_plans()[pme_unit];
  auto& st = detail_::pme_objs()[pme_unit];

#  if defined(TINKER_GPU_SINGLE)
  cufftExecC2C(iplan, reinterpret_cast<cufftComplex*>(st.qgrid),
               reinterpret_cast<cufftComplex*>(st.qgrid), CUFFT_INVERSE);
#  elif defined(TINKER_GPU_DOUBLE)
  cufftExecZ2Z(iplan, reinterpret_cast<cufftDoubleComplex*>(st.qgrid),
               reinterpret_cast<cufftDoubleComplex*>(st.qgrid), CUFFT_INVERSE);
#  else
  static_assert(false, "");
#  endif
}
}
TINKER_NAMESPACE_END

#endif
