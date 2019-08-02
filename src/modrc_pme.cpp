#include "mod_md.h"
#include "mod_pme.h"
#include "util_array.h"
#include "util_math.h"
#include "util_potent.h"
#include "util_switch.h"
#include <ext/tinker/tinker_mod.h>
#include <ext/tinker/tinker_rt.h>

TINKER_NAMESPACE_BEGIN
pme_t& pme_obj(PMEUnit pme_u) {
#if TINKER_DEBUG
  return PMEUnit::all_objs().at(pme_u.unit());
#else
  return PMEUnit::all_objs()[pme_u.unit()];
#endif
}

pme_t* pme_deviceptr(PMEUnit pme_u) {
  int u = pme_u.unit();
#if TINKER_DEBUG
  return PMEUnit::all_deviceptrs().at(pme_u.unit());
#else
  return PMEUnit::all_deviceptrs()[pme_u.unit()];
#endif
}

static void pme_op_dealloc_(int pu) {
  pme_t& st = pme_obj(pu);
  pme_t* dptr = pme_deviceptr(pu);

  check_rt(cudaFree(st.igrid));
  check_rt(cudaFree(st.bsmod1));
  check_rt(cudaFree(st.bsmod2));
  check_rt(cudaFree(st.bsmod3));
  check_rt(cudaFree(st.qgrid));

  check_rt(cudaFree(dptr));
}

static void pme_op_alloc_(PMEUnit& unit, double aewald, int nfft1, int nfft2,
                          int nfft3, int bsorder, bool unique) {
  int count = 0;
  int first = -1;
  const double eps = 1.0e-6;
  int idx;
  for (idx = 0; idx < PMEUnit::all_objs().size(); ++idx) {
    auto& st = PMEUnit::all_objs()[idx];
    if (std::abs(aewald - st.aewald) < eps && nfft1 == st.nfft1 &&
        nfft2 == st.nfft2 && nfft3 == st.nfft3 && bsorder == st.bsorder) {
      ++count;
      if (count == 1)
        first = idx;
    }
  }

  if (count == 0 || unique == true) {
    PMEUnit::all_objs().emplace_back(pme_t());
    PMEUnit::all_deviceptrs().emplace_back(nullptr);
    auto& st = PMEUnit::all_objs()[idx];
    auto& dptr = PMEUnit::all_deviceptrs()[idx];

    const size_t rs = sizeof(real);
    size_t size;

    size = 3 * n * sizeof(int);
    check_rt(cudaMalloc(&st.igrid, size));
    // see also subroutine moduli in pmestuf.f
    check_rt(cudaMalloc(&st.bsmod1, rs * nfft1));
    check_rt(cudaMalloc(&st.bsmod2, rs * nfft2));
    check_rt(cudaMalloc(&st.bsmod3, rs * nfft3));
    size = nfft1 * nfft2 * nfft3 * rs;
    check_rt(cudaMalloc(&st.qgrid, 2 * size));

    size = sizeof(pme_t);
    check_rt(cudaMalloc(&dptr, size));

    st.aewald = aewald;
    st.nfft1 = nfft1;
    st.nfft2 = nfft2;
    st.nfft3 = nfft3;
    st.bsorder = bsorder;
  } else {
    idx = first;
  }

  assert(PMEUnit::all_objs().size() == PMEUnit::all_deviceptrs().size());
  unit = idx;
}

static void pme_op_copyin_(PMEUnit unit) {
  if (unit < 0)
    return;

  pme_t& st = pme_obj(unit);
  pme_t* dptr = pme_deviceptr(unit);

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

  size_t size = sizeof(pme_t);
  check_rt(cudaMemcpy(dptr, &st, size, cudaMemcpyHostToDevice));
}
TINKER_NAMESPACE_END

#include "gpu/e_mpole.h"
#include "gpu/e_polar.h"

TINKER_NAMESPACE_BEGIN

int use_ewald() { return limits::use_ewald; }

void pme_init(int vers) {
  if (!use_ewald())
    return;

  rpole_to_cmp();

  if (vir_m) {
    zero_array(vir_m, 9);
  }
}

void pme_data(rc_op op) {
  if (!use_ewald())
    return;

  if (op & rc_dealloc) {
    assert(PMEUnit::all_objs().size() == PMEUnit::all_deviceptrs().size());
    int idx = 0;
    while (idx < PMEUnit::all_objs().size()) {
      pme_op_dealloc_(idx);
      ++idx;
    }
    PMEUnit::all_objs().clear();
    PMEUnit::all_deviceptrs().clear();

    check_rt(cudaFree(cmp));
    check_rt(cudaFree(fmp));
    check_rt(cudaFree(cphi));
    check_rt(cudaFree(fphi));

    if (use_potent(polar_term)) {
      check_rt(cudaFree(fuind));
      check_rt(cudaFree(fuinp));
      check_rt(cudaFree(fdip_phi1));
      check_rt(cudaFree(fdip_phi2));
      check_rt(cudaFree(cphidp));
      check_rt(cudaFree(fphidp));
    }

    check_rt(cudaFree(vir_m));
  }

  if (op & rc_alloc) {
    assert(PMEUnit::all_objs().size() == 0);
    assert(PMEUnit::all_deviceptrs().size() == 0);

    const size_t rs = sizeof(real);

    check_rt(cudaMalloc(&cmp, 10 * n * rs));
    check_rt(cudaMalloc(&fmp, 10 * n * rs));
    check_rt(cudaMalloc(&cphi, 10 * n * rs));
    check_rt(cudaMalloc(&fphi, 20 * n * rs));

    if (use_potent(polar_term)) {
      check_rt(cudaMalloc(&fuind, 3 * n * rs));
      check_rt(cudaMalloc(&fuinp, 3 * n * rs));
      check_rt(cudaMalloc(&fdip_phi1, 10 * n * rs));
      check_rt(cudaMalloc(&fdip_phi2, 10 * n * rs));
      check_rt(cudaMalloc(&cphidp, 10 * n * rs));
      check_rt(cudaMalloc(&fphidp, 20 * n * rs));

      // if (vir_m), it implies use virial and use epolar
      if (use_data & calc::virial)
        check_rt(cudaMalloc(&vir_m, 9 * rs));
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
      if (use_data & calc::virial) {
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

  if (op & rc_init) {
    switch_cut_off(switch_ewald, ewald_switch_cut, ewald_switch_off);

    pme_op_copyin_(epme_unit);
    pme_op_copyin_(ppme_unit);
    pme_op_copyin_(pvpme_unit);
    pme_op_copyin_(dpme_unit);
  }

  fft_data(op);
}
TINKER_NAMESPACE_END
