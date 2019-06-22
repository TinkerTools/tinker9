#ifdef TINKER_HOSTONLY

#  include "gpu/decl_dataop.h"
#  include "gpu/decl_pme.h"
#  include "rc.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
std::vector<fft_plan_t>& fft_plans();
}
TINKER_NAMESPACE_END

//======================================================================

TINKER_NAMESPACE_BEGIN
namespace gpu {
void fft_data(int op) {
  if (op & op_dealloc) {
    int idx = 0;
    while (idx < fft_plans().size()) {
      auto& ps = fft_plans()[idx];

#  if defined(TINKER_GPU_SINGLE)
      fftwf_destroy_plan(ps.planf);
      fftwf_destroy_plan(ps.planb);
#  elif defined(TINKER_GPU_DOUBLE)
      fftw_destroy_plan(ps.planf);
      fftw_destroy_plan(ps.planb);
#  else
      static_assert(false, "");
#  endif

      ++idx;
    }

    fft_plans().clear();
  }

  if (op & op_alloc) {
    assert(fft_plans().size() == 0);

    const size_t size = detail_::pme_objs().size();
    fft_plans().resize(size, fft_plan_t());
  }

  if (op & op_copyin) {
    int idx = 0;
    for (fft_plan_t& iplan : fft_plans()) {
      auto& st = pme_obj(idx);
      const int nfft1 = st.nfft1;
      const int nfft2 = st.nfft2;
      const int nfft3 = st.nfft3;
      const int ifront = -1;
      const int iback = 1;
      const unsigned int iguess = 0;

#  if defined(TINKER_GPU_SINGLE)
      iplan.planf = fftwf_plan_dft_3d(
          nfft1, nfft2, nfft3, reinterpret_cast<fftwf_complex*>(st.qgrid),
          reinterpret_cast<fftwf_complex*>(st.qgrid), ifront, iguess);
      iplan.planb = fftwf_plan_dft_3d(
          nfft1, nfft2, nfft3, reinterpret_cast<fftwf_complex*>(st.qgrid),
          reinterpret_cast<fftwf_complex*>(st.qgrid), iback, iguess);
#  elif defined(TINKER_GPU_DOUBLE)
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

void fftfront(int pme_unit) {
  fft_plan_t iplan = fft_plans()[pme_unit];
  auto& st = pme_obj(pme_unit);

#  if defined(TINKER_GPU_SINGLE)
  fftwf_execute_dft(iplan.planf, reinterpret_cast<fftwf_complex*>(st.qgrid),
                    reinterpret_cast<fftwf_complex*>(st.qgrid));
#  elif defined(TINKER_GPU_DOUBLE)
  fftw_execute_dft(iplan.planf, reinterpret_cast<fftw_complex*>(st.qgrid),
                   reinterpret_cast<fftw_complex*>(st.qgrid));
#  else
  static_assert(false, "");
#  endif
}

void fftback(int pme_unit) {
  fft_plan_t iplan = fft_plans()[pme_unit];
  auto& st = pme_obj(pme_unit);

#  if defined(TINKER_GPU_SINGLE)
  fftwf_execute_dft(iplan.planb, reinterpret_cast<fftwf_complex*>(st.qgrid),
                    reinterpret_cast<fftwf_complex*>(st.qgrid));
#  elif defined(TINKER_GPU_DOUBLE)
  fftw_execute_dft(iplan.planb, reinterpret_cast<fftw_complex*>(st.qgrid),
                   reinterpret_cast<fftw_complex*>(st.qgrid));
#  else
  static_assert(false, "");
#  endif
}
}
TINKER_NAMESPACE_END

#endif
