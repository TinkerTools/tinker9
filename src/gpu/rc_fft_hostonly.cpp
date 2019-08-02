#ifdef TINKER_HOSTONLY

#  include "mod_pme.h"
#  include "util_rt.h"

TINKER_NAMESPACE_BEGIN
extern std::vector<FFTPlan>& fft_plans();

void fft_data(rc_op op) {
  if (op & rc_dealloc) {
    int idx = 0;
    while (idx < fft_plans().size()) {
      auto& ps = fft_plans()[idx];

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

    fft_plans().clear();
  }

  if (op & rc_alloc) {
    assert(fft_plans().size() == 0);

    const size_t size = PMEUnit::all_objs().size();
    fft_plans().resize(size, FFTPlan());
  }

  if (op & rc_init) {
    int idx = 0;
    for (FFTPlan& iplan : fft_plans()) {
      auto& st = pme_obj(idx);
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
  FFTPlan iplan = fft_plans()[pme_u.unit()];
  auto& st = pme_obj(pme_u.unit());

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
  FFTPlan iplan = fft_plans()[pme_u.unit()];
  auto& st = pme_obj(pme_u.unit());

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
