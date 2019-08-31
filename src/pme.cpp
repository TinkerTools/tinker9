#include "pme.h"
#include "array.h"
#include "elec.h"
#include "ext/tinker/detail/ewald.hh"
#include "ext/tinker/detail/pme.hh"
#include "mathfunc.h"
#include "md.h"
#include "potent.h"
#include "switch.h"
#include "tinker_rt.h"

TINKER_NAMESPACE_BEGIN
bool PME::Params::operator==(const Params& st) const {
  const double eps = 1.0e-6;
  bool ans = std::fabs(aewald - st.aewald) < eps && nfft1 == st.nfft1 &&
      nfft2 == st.nfft2 && nfft3 == st.nfft3 && bsorder == st.bsorder;
  return ans;
}

PME::Params::Params(real a, int n1, int n2, int n3, int o)
    : aewald(a)
    , nfft1(n1)
    , nfft2(n2)
    , nfft3(n3)
    , bsorder(o) {}

void PME::set_params(const PME::Params& p) {
  aewald = p.aewald;
  nfft1 = p.nfft1;
  nfft2 = p.nfft2;
  nfft3 = p.nfft3;
  bsorder = p.bsorder;
}

PME::Params PME::get_params() const {
  Params p0(aewald, nfft1, nfft2, nfft3, bsorder);
  return p0;
}

bool PME::operator==(const Params& p) const { return get_params() == p; }

PME::~PME() {
  dealloc_bytes(bsmod1);
  dealloc_bytes(bsmod2);
  dealloc_bytes(bsmod3);
  dealloc_bytes(qgrid);
}

static void pme_op_alloc_(PMEUnit& unit, const PME::Params& p, bool unique) {
  unit.close();
  for (PMEUnit idx = 0; idx < PMEUnit::size(); idx = idx + 1) {
    if (*idx == p)
      unit = idx;
  }

  if (!unit.valid() || unique == true) {
    unit = PMEUnit::open();
    auto& st = *unit;
    const size_t rs = sizeof(real);
    size_t size;

    // see also subroutine moduli in pmestuf.f
    alloc_bytes(&st.bsmod1, rs * p.nfft1);
    alloc_bytes(&st.bsmod2, rs * p.nfft2);
    alloc_bytes(&st.bsmod3, rs * p.nfft3);
    size = p.nfft1 * p.nfft2 * p.nfft3 * rs;
    alloc_bytes(&st.qgrid, 2 * size);

    st.set_params(p);
  }
}

static void pme_op_copyin_(PMEUnit unit) {
  if (!unit.valid())
    return;

  auto& st = *unit;

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

  unit.init_deviceptr(st);
}
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
void pme_init(int vers) {
  if (!use_ewald())
    return;

  rpole_to_cmp();

  if (vir_m_handle.valid())
    vir_m_handle->zero();
}

static void pme_data1_(rc_op op) {
  if (op & rc_dealloc) {
    PMEUnit::clear();

    if (use_potent(polar_term)) {
      vir_m_handle.close();
    }

    epme_unit.close();
    ppme_unit.close();
    pvpme_unit.close();
    dpme_unit.close();
  }

  if (op & rc_alloc) {
    assert(PMEUnit::size() == 0);

    cmp_vec.reserve(10 * n);
    fmp_vec.reserve(10 * n);
    cphi_vec.reserve(10 * n);
    fphi_vec.reserve(20 * n);

    if (use_potent(polar_term)) {
      fuind_vec.reserve(3 * n);
      fuinp_vec.reserve(3 * n);
      fdip_phi1_vec.reserve(10 * n);
      fdip_phi2_vec.reserve(10 * n);
      cphidp_vec.reserve(10 * n);
      fphidp_vec.reserve(20 * n);

      if (rc_flag & calc::virial) {
        vir_m_handle = Virial::open();
        vir_m_handle->alloc(n);
      }
    }

    bool unique_grids = false;

    // electrostatics
    epme_unit.close();
    if (use_potent(mpole_term)) {
      unique_grids = false;
      PME::Params p(ewald::aeewald, pme::nefft1, pme::nefft2, pme::nefft3,
                    pme::bseorder);
      pme_op_alloc_(epme_unit, p, unique_grids);
    }

    // polarization
    ppme_unit.close();
    pvpme_unit.close();
    if (use_potent(polar_term)) {
      PME::Params p(ewald::apewald, pme::nefft1, pme::nefft2, pme::nefft3,
                    pme::bsporder);
      pme_op_alloc_(ppme_unit, p, unique_grids);
      if (rc_flag & calc::virial) {
        unique_grids = true;
        pme_op_alloc_(pvpme_unit, p, unique_grids);
      }
    }

    // dispersion
    dpme_unit.close();
    if (false) {
      unique_grids = false;
      PME::Params p(ewald::adewald, pme::ndfft1, pme::ndfft2, pme::ndfft3,
                    pme::bsdorder);
      pme_op_alloc_(dpme_unit, p, unique_grids);
    }
  }

  if (op & rc_init) {
    switch_cut_off(switch_ewald, ewald_switch_cut, ewald_switch_off);

    pme_op_copyin_(epme_unit);
    pme_op_copyin_(ppme_unit);
    pme_op_copyin_(pvpme_unit);
    pme_op_copyin_(dpme_unit);
  }
}

void pme_data(rc_op op) {
  if (!use_ewald())
    return;

  rc_man pme42_{pme_data1_, op};
  rc_man fft42_{fft_data, op};
}
TINKER_NAMESPACE_END
