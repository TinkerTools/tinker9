#ifndef TINKER_HOSTONLY

#  include "mod_pme.h"
#  include "util_rt.h"

TINKER_NAMESPACE_BEGIN
void fft_data(rc_op op) {
  if (op & rc_dealloc) {
    int idx = 0;
    while (idx < FFTPlanUnit::size()) {
      FFTPlanUnit u = idx;
      cufftDestroy(u.obj());
      ++idx;
    }
    FFTPlanUnit::clear();
  }

  if (op & rc_alloc) {
    assert(FFTPlanUnit::size() == 0);

    const size_t size = PMEUnit::size();
    FFTPlanUnit::resize(size);
  }

  if (op & rc_init) {
#  if defined(TINKER_SINGLE_PRECISION)
    const cufftType typ = CUFFT_C2C;
#  elif defined(TINKER_DOUBLE_PRECISION)
    const cufftType typ = CUFFT_Z2Z;
#  else
    static_assert(false, "");
#  endif

    int idx = 0;
    while (idx < FFTPlanUnit::size()) {
      FFTPlanUnit plan_u = idx;
      PMEUnit pme_u = idx;
      auto& iplan = plan_u.obj();
      auto& st = pme_u.obj();

      check_rt(cufftPlan3d(&iplan, st.nfft1, st.nfft2, st.nfft3, typ));
      ++idx;
    }
  }
}

void fftfront(PMEUnit pme_u) {
  FFTPlanUnit iplan_u = static_cast<int>(pme_u);
  auto& iplan = iplan_u.obj();
  auto& st = pme_u.obj();

#  if defined(TINKER_SINGLE_PRECISION)
  cufftExecC2C(iplan, reinterpret_cast<cufftComplex*>(st.qgrid),
               reinterpret_cast<cufftComplex*>(st.qgrid), CUFFT_FORWARD);
#  elif defined(TINKER_DOUBLE_PRECISION)
  cufftExecZ2Z(iplan, reinterpret_cast<cufftDoubleComplex*>(st.qgrid),
               reinterpret_cast<cufftDoubleComplex*>(st.qgrid), CUFFT_FORWARD);
#  else
  static_assert(false, "");
#  endif
}

void fftback(PMEUnit pme_u) {
  FFTPlanUnit iplan_u = static_cast<int>(pme_u);
  auto& iplan = iplan_u.obj();
  auto& st = pme_u.obj();

#  if defined(TINKER_SINGLE_PRECISION)
  cufftExecC2C(iplan, reinterpret_cast<cufftComplex*>(st.qgrid),
               reinterpret_cast<cufftComplex*>(st.qgrid), CUFFT_INVERSE);
#  elif defined(TINKER_DOUBLE_PRECISION)
  cufftExecZ2Z(iplan, reinterpret_cast<cufftDoubleComplex*>(st.qgrid),
               reinterpret_cast<cufftDoubleComplex*>(st.qgrid), CUFFT_INVERSE);
#  else
  static_assert(false, "");
#  endif
}
TINKER_NAMESPACE_END

#endif
