#include "array.h"
#include "mathfunc.h"
#include "md.h"
#include "pme.h"
#include "switch.h"
#include "util_potent.h"
#include <ext/tinker/tinker_mod.h>
#include <ext/tinker/tinker_rt.h>

TINKER_NAMESPACE_BEGIN
static void pme_op_dealloc_(PMEUnit pu) {
  auto& st = pu.obj();
  auto* dptr = pu.deviceptr();

  dealloc_bytes(st.igrid);
  dealloc_bytes(st.bsmod1);
  dealloc_bytes(st.bsmod2);
  dealloc_bytes(st.bsmod3);
  dealloc_bytes(st.qgrid);

  dealloc_bytes(dptr);
}

static void pme_op_alloc_(PMEUnit& unit, double aewald, int nfft1, int nfft2,
                          int nfft3, int bsorder, bool unique) {
  int count = 0;
  int first = -1;
  const double eps = 1.0e-6;
  int idx;
  PMEUnit st_u;
  for (idx = 0; idx < PMEUnit::size(); ++idx) {
    st_u = idx;
    auto& st = st_u.obj();
    if (std::abs(aewald - st.aewald) < eps && nfft1 == st.nfft1 &&
        nfft2 == st.nfft2 && nfft3 == st.nfft3 && bsorder == st.bsorder) {
      ++count;
      if (count == 1)
        first = idx;
    }
  }

  if (count == 0 || unique == true) {
    PMEUnit::emplace_back(PME());
    st_u = idx;
    auto& st = st_u.obj();
    PME*& dptr = st_u.deviceptr();

    const size_t rs = sizeof(real);
    size_t size;

    size = 3 * n * sizeof(int);
    alloc_bytes(&st.igrid, size);
    // see also subroutine moduli in pmestuf.f
    alloc_bytes(&st.bsmod1, rs * nfft1);
    alloc_bytes(&st.bsmod2, rs * nfft2);
    alloc_bytes(&st.bsmod3, rs * nfft3);
    size = nfft1 * nfft2 * nfft3 * rs;
    alloc_bytes(&st.qgrid, 2 * size);

    size = sizeof(PME);
    alloc_bytes(&dptr, size);

    st.aewald = aewald;
    st.nfft1 = nfft1;
    st.nfft2 = nfft2;
    st.nfft3 = nfft3;
    st.bsorder = bsorder;
  } else {
    idx = first;
  }

  unit = idx;
}

static void pme_op_copyin_(PMEUnit unit) {
  if (unit < 0)
    return;

  auto& st = unit.obj();
  auto* dptr = unit.deviceptr();

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

  size_t size = sizeof(PME);
  copyin_bytes(dptr, &st, size);
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

    int idx = 0;
    while (idx < PMEUnit::size()) {
      PMEUnit pu = idx;
      pme_op_dealloc_(pu);
      ++idx;
    }
    PMEUnit::clear();

    dealloc_bytes(cmp);
    dealloc_bytes(fmp);
    dealloc_bytes(cphi);
    dealloc_bytes(fphi);

    if (use_potent(polar_term)) {
      dealloc_bytes(fuind);
      dealloc_bytes(fuinp);
      dealloc_bytes(fdip_phi1);
      dealloc_bytes(fdip_phi2);
      dealloc_bytes(cphidp);
      dealloc_bytes(fphidp);
    }

    dealloc_bytes(vir_m);
  }

  if (op & rc_alloc) {
    assert(PMEUnit::size() == 0);

    const size_t rs = sizeof(real);

    alloc_bytes(&cmp, 10 * n * rs);
    alloc_bytes(&fmp, 10 * n * rs);
    alloc_bytes(&cphi, 10 * n * rs);
    alloc_bytes(&fphi, 20 * n * rs);

    if (use_potent(polar_term)) {
      alloc_bytes(&fuind, 3 * n * rs);
      alloc_bytes(&fuinp, 3 * n * rs);
      alloc_bytes(&fdip_phi1, 10 * n * rs);
      alloc_bytes(&fdip_phi2, 10 * n * rs);
      alloc_bytes(&cphidp, 10 * n * rs);
      alloc_bytes(&fphidp, 20 * n * rs);

      // if (vir_m), it implies use virial and use epolar
      if (rc_flag & calc::virial)
        alloc_bytes(&vir_m, 9 * rs);
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
      if (rc_flag & calc::virial) {
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
