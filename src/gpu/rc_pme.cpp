#include "gpu/decl_mdstate.h"
#include "gpu/decl_pme.h"
#include "gpu/decl_potent.h"
#include "gpu/decl_switch.h"
#include "gpu/rc.h"
#include "util/math.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
namespace detail_ {
std::vector<pme_st>& pme_objs() {
  static std::vector<pme_st> objs;
  return objs;
}

std::vector<pme_st*>& pme_deviceptrs() {
  static std::vector<pme_st*> ptrs;
  return ptrs;
}
}

pme_st& pme_obj(int pme_unit) {
#if TINKER_DEBUG
  return detail_::pme_objs().at(pme_unit);
#else
  return detail_::pme_objs()[pme_unit];
#endif
}

pme_st* pme_deviceptr(int pme_unit) {
#if TINKER_DEBUG
  return detail_::pme_deviceptrs().at(pme_unit);
#else
  return detail_::pme_deviceptrs()[pme_unit];
#endif
}

static void pme_op_dealloc_(int pu) {
  pme_st& st = pme_obj(pu);
  pme_st* dptr = pme_deviceptr(pu);

  check_cudart(cudaFree(st.igrid));
  check_cudart(cudaFree(st.bsmod1));
  check_cudart(cudaFree(st.bsmod2));
  check_cudart(cudaFree(st.bsmod3));
  check_cudart(cudaFree(st.qgrid));

  check_cudart(cudaFree(dptr));
}

static int pme_op_alloc_(int& unit, double aewald, int nfft1, int nfft2,
                         int nfft3, int bsorder, bool unique) {
  int count = 0;
  int first = -1;
  const double eps = 1.0e-6;
  int idx;
  for (idx = 0; idx < detail_::pme_objs().size(); ++idx) {
    auto& st = detail_::pme_objs()[idx];
    if (std::abs(aewald - st.aewald) < eps && nfft1 == st.nfft1 &&
        nfft2 == st.nfft2 && nfft3 == st.nfft3 && bsorder == st.bsorder) {
      ++count;
      if (count == 1)
        first = idx;
    }
  }

  if (count == 0 || unique == true) {
    detail_::pme_objs().emplace_back(pme_st());
    detail_::pme_deviceptrs().emplace_back(nullptr);
    auto& st = detail_::pme_objs()[idx];
    auto& dptr = detail_::pme_deviceptrs()[idx];

    const size_t rs = sizeof(real);
    size_t size;

    size = 3 * n * sizeof(int);
    check_cudart(cudaMalloc(&st.igrid, size));
    // see also subroutine moduli in pmestuf.f
    check_cudart(cudaMalloc(&st.bsmod1, rs * nfft1));
    check_cudart(cudaMalloc(&st.bsmod2, rs * nfft2));
    check_cudart(cudaMalloc(&st.bsmod3, rs * nfft3));
    size = nfft1 * nfft2 * nfft3 * rs;
    check_cudart(cudaMalloc(&st.qgrid, 2 * size));

    size = sizeof(pme_st);
    check_cudart(cudaMalloc(&dptr, size));

    st.aewald = aewald;
    st.nfft1 = nfft1;
    st.nfft2 = nfft2;
    st.nfft3 = nfft3;
    st.bsorder = bsorder;
  } else {
    idx = first;
  }

  assert(detail_::pme_objs().size() == detail_::pme_deviceptrs().size());
  unit = idx;
}

static void pme_op_copyin_(int unit) {
  if (unit < 0)
    return;

  pme_st& st = pme_obj(unit);
  pme_st* dptr = pme_deviceptr(unit);

  // This code assumes that the FFT grids of an energy term will not change in a
  // calculation.
  int maxfft = max_of(st.nfft1, st.nfft2, st.nfft3);
  std::vector<double> array(st.bsorder);
  std::vector<double> bsarray(maxfft);
  double x = 0;
  TINKER_RT(bspline)(&x, &st.bsorder, array.data());
  for (int i = 0; i < maxfft; ++i) {
    bsarray[i] = 0;
  }
  assert(st.bsorder + 1 <= maxfft);
  for (int i = 0; i < st.bsorder; ++i) {
    bsarray[i + 1] = array[i];
  }
  std::vector<double> bsmodbuf(maxfft);
  TINKER_RT(dftmod)(bsmodbuf.data(), bsarray.data(), &st.nfft1, &st.bsorder);
  copyin_array(st.bsmod1, bsmodbuf.data(), st.nfft1);
  TINKER_RT(dftmod)(bsmodbuf.data(), bsarray.data(), &st.nfft2, &st.bsorder);
  copyin_array(st.bsmod2, bsmodbuf.data(), st.nfft2);
  TINKER_RT(dftmod)(bsmodbuf.data(), bsarray.data(), &st.nfft3, &st.bsorder);
  copyin_array(st.bsmod3, bsmodbuf.data(), st.nfft3);

  size_t size = sizeof(pme_st);
  check_cudart(cudaMemcpy(dptr, &st, size, cudaMemcpyHostToDevice));
}

int epme_unit; // electrostatic
int ppme_unit; // polarization
int dpme_unit; // dispersion

int pvpme_unit; // polarization virial

double ewald_switch_cut, ewald_switch_off;
}
TINKER_NAMESPACE_END

#include "gpu/e_mpole.h"
#include "gpu/e_polar.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
real (*cmp)[10];
real (*fmp)[10];
real (*cphi)[10];
real (*fphi)[20];

real (*fuind)[3];
real (*fuinp)[3];
real (*fdip_phi1)[10];
real (*fdip_phi2)[10];
real (*cphidp)[10];
real (*fphidp)[20];

real* vir_m;

int use_ewald() { return limits::use_ewald; }

void pme_init(int vers) {
  if (!use_ewald())
    return;

  rpole_to_cmp();

  if (vir_m) {
    zero_array(vir_m, 9);
  }
}

void pme_data(rc_t rc) {
  if (!use_ewald())
    return;

  if (rc & rc_dealloc) {
    assert(detail_::pme_objs().size() == detail_::pme_deviceptrs().size());
    int idx = 0;
    while (idx < detail_::pme_objs().size()) {
      pme_op_dealloc_(idx);
      ++idx;
    }
    detail_::pme_objs().clear();
    detail_::pme_deviceptrs().clear();

    check_cudart(cudaFree(cmp));
    check_cudart(cudaFree(fmp));
    check_cudart(cudaFree(cphi));
    check_cudart(cudaFree(fphi));

    if (use_potent(polar_term)) {
      check_cudart(cudaFree(fuind));
      check_cudart(cudaFree(fuinp));
      check_cudart(cudaFree(fdip_phi1));
      check_cudart(cudaFree(fdip_phi2));
      check_cudart(cudaFree(cphidp));
      check_cudart(cudaFree(fphidp));
    }

    check_cudart(cudaFree(vir_m));
  }

  if (rc & rc_alloc) {
    assert(detail_::pme_objs().size() == 0);
    assert(detail_::pme_deviceptrs().size() == 0);

    const size_t rs = sizeof(real);

    check_cudart(cudaMalloc(&cmp, 10 * n * rs));
    check_cudart(cudaMalloc(&fmp, 10 * n * rs));
    check_cudart(cudaMalloc(&cphi, 10 * n * rs));
    check_cudart(cudaMalloc(&fphi, 20 * n * rs));

    if (use_potent(polar_term)) {
      check_cudart(cudaMalloc(&fuind, 3 * n * rs));
      check_cudart(cudaMalloc(&fuinp, 3 * n * rs));
      check_cudart(cudaMalloc(&fdip_phi1, 10 * n * rs));
      check_cudart(cudaMalloc(&fdip_phi2, 10 * n * rs));
      check_cudart(cudaMalloc(&cphidp, 10 * n * rs));
      check_cudart(cudaMalloc(&fphidp, 20 * n * rs));

      // if (vir_m), it implies use virial and use epolar
      if (use_data & use_virial)
        check_cudart(cudaMalloc(&vir_m, 9 * rs));
      else
        vir_m = nullptr;
    }

    bool unique_grids = false;

    // electrostatics
    epme_unit = -1;
    if (use_potent(mpole_term)) {
      unique_grids = false;
      pme_op_alloc_(epme_unit, ewald::aeewald, pme::nefft1, pme::nefft2,
                    pme::nefft3, pme::bseorder, unique_grids);
    }

    // polarization
    ppme_unit = -1;
    pvpme_unit = -1;
    if (use_potent(polar_term)) {
      pme_op_alloc_(ppme_unit, ewald::apewald, pme::nefft1, pme::nefft2,
                    pme::nefft3, pme::bsporder, unique_grids);
      if (use_data & use_virial) {
        unique_grids = true;
        pme_op_alloc_(pvpme_unit, ewald::apewald, pme::nefft1, pme::nefft2,
                      pme::nefft3, pme::bsporder, unique_grids);
      }
    }

    // dispersion
    dpme_unit = -1;
    if (false) {
      unique_grids = false;
      pme_op_alloc_(dpme_unit, ewald::adewald, pme::ndfft1, pme::ndfft2,
                    pme::ndfft3, pme::bsdorder, unique_grids);
    }
  }

  if (rc & rc_copyin) {
    switch_cut_off(switch_ewald, ewald_switch_cut, ewald_switch_off);

    pme_op_copyin_(epme_unit);
    pme_op_copyin_(ppme_unit);
    pme_op_copyin_(pvpme_unit);
    pme_op_copyin_(dpme_unit);
  }

  fft_data(rc);
}
}
TINKER_NAMESPACE_END
