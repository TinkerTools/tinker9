#include "pme.h"
#include "rt.h"

#ifdef TINKER_HOST
TINKER_NAMESPACE_BEGIN
void fft_data(rc_op op) {
  if (op & rc_dealloc) {
    int idx = 0;
    while (idx < FFTPlanUnit::size()) {
      FFTPlanUnit u = idx;
      auto& ps = u.obj();
#  if defined(TINKER_SINGLE_PRECISION)
      fftwf_destroy_plan(ps.planf);
      fftwf_destroy_plan(ps.planb);
#  elif defined(TINKER_DOUBLE_PRECISION)
      fftw_destroy_plan(ps.planf);
      fftw_destroy_plan(ps.planb);
#  else
      static_assert(false, "");
#  endif
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
    int idx = 0;
    while (idx < FFTPlanUnit::size()) {
      FFTPlanUnit plan_u = idx;
      PMEUnit pme_u = idx;
      auto& iplan = plan_u.obj();
      auto& st = pme_u.obj();

      const int nfft1 = st.nfft1;
      const int nfft2 = st.nfft2;
      const int nfft3 = st.nfft3;
      const int ifront = -1;
      const int iback = 1;
      const unsigned int iguess = 0;

#  if defined(TINKER_SINGLE_PRECISION)
      iplan.planf = fftwf_plan_dft_3d(
          nfft1, nfft2, nfft3, reinterpret_cast<fftwf_complex*>(st.qgrid),
          reinterpret_cast<fftwf_complex*>(st.qgrid), ifront, iguess);
      iplan.planb = fftwf_plan_dft_3d(
          nfft1, nfft2, nfft3, reinterpret_cast<fftwf_complex*>(st.qgrid),
          reinterpret_cast<fftwf_complex*>(st.qgrid), iback, iguess);
#  elif defined(TINKER_DOUBLE_PRECISION)
      iplan.planf = fftw_plan_dft_3d(
          nfft1, nfft2, nfft3, reinterpret_cast<fftw_complex*>(st.qgrid),
          reinterpret_cast<fftw_complex*>(st.qgrid), ifront, iguess);
      iplan.planb = fftw_plan_dft_3d(
          nfft1, nfft2, nfft3, reinterpret_cast<fftw_complex*>(st.qgrid),
          reinterpret_cast<fftw_complex*>(st.qgrid), iback, iguess);
#  else
      static_assert(false, "");
#  endif
      ++idx;
    }
  }
}

void fftfront(PMEUnit pme_u) {
  FFTPlanUnit iplan_u = static_cast<int>(pme_u);
  auto& iplan = iplan_u.obj();
  auto& st = pme_u.obj();

#  if defined(TINKER_SINGLE_PRECISION)
  fftwf_execute_dft(iplan.planf, reinterpret_cast<fftwf_complex*>(st.qgrid),
                    reinterpret_cast<fftwf_complex*>(st.qgrid));
#  elif defined(TINKER_DOUBLE_PRECISION)
  fftw_execute_dft(iplan.planf, reinterpret_cast<fftw_complex*>(st.qgrid),
                   reinterpret_cast<fftw_complex*>(st.qgrid));
#  else
  static_assert(false, "");
#  endif
}

void fftback(PMEUnit pme_u) {
  FFTPlanUnit iplan_u = static_cast<int>(pme_u);
  auto& iplan = iplan_u.obj();
  auto& st = pme_u.obj();

#  if defined(TINKER_SINGLE_PRECISION)
  fftwf_execute_dft(iplan.planb, reinterpret_cast<fftwf_complex*>(st.qgrid),
                    reinterpret_cast<fftwf_complex*>(st.qgrid));
#  elif defined(TINKER_DOUBLE_PRECISION)
  fftw_execute_dft(iplan.planb, reinterpret_cast<fftw_complex*>(st.qgrid),
                   reinterpret_cast<fftw_complex*>(st.qgrid));
#  else
  static_assert(false, "");
#  endif
}
TINKER_NAMESPACE_END

#endif
