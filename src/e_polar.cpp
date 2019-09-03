#include "e_polar.h"

#include "ext/tinker/detail/inform.hh"
#include "ext/tinker/detail/polar.hh"
#include "ext/tinker/detail/polpot.hh"
#include "ext/tinker/detail/units.hh"
#include "io_fort_str.h"
#include "io_print.h"
#include "md.h"
#include "pme.h"
#include "potent.h"
#include "switch.h"

TINKER_NAMESPACE_BEGIN
void epolar_data(rc_op op) {
  if (!use_potent(polar_term))
    return;

  if (op & rc_dealloc) {
    device_array::deallocate(polarity, thole, pdamp, polarity_inv);

    ep_handle.dealloc();

    device_array::deallocate(ufld, dufld);
    device_array::deallocate(work01_, work02_, work03_, work04_, work05_,
                             work06_, work07_, work08_, work09_, work10_);
  }

  if (op & rc_alloc) {
    device_array::allocate(n, &polarity, &thole, &pdamp, &polarity_inv);

    ep_handle.alloc(n);

    if (rc_flag & calc::grad) {
      device_array::allocate(n, &ufld, &dufld);
    } else {
      ufld = nullptr;
      dufld = nullptr;
    }

    device_array::allocate(n, &work01_, &work02_, &work03_, &work04_, &work05_,
                           &work06_, &work07_, &work08_, &work09_, &work10_);
  }

  if (op & rc_init) {
    if (use_ewald()) {
      epolar_electyp = elec_t::ewald;
    } else {
      epolar_electyp = elec_t::coulomb;
    }

    if (epolar_electyp == elec_t::coulomb)
      switch_cut_off(switch_mpole, mpole_switch_cut, mpole_switch_off);

    u1scale = polpot::u1scale;
    u2scale = polpot::u2scale;
    u3scale = polpot::u3scale;
    u4scale = polpot::u4scale;

    d1scale = polpot::d1scale;
    d2scale = polpot::d2scale;
    d3scale = polpot::d3scale;
    d4scale = polpot::d4scale;

    p2scale = polpot::p2scale;
    p3scale = polpot::p3scale;
    p4scale = polpot::p4scale;
    p5scale = polpot::p5scale;

    p2iscale = polpot::p2iscale;
    p3iscale = polpot::p3iscale;
    p4iscale = polpot::p4iscale;
    p5iscale = polpot::p5iscale;

    udiag = polpot::udiag;

    // see also polmin in induce.f
    const double polmin = 0.00000001;
    std::vector<double> pinvbuf(n);
    for (int i = 0; i < n; ++i) {
      pinvbuf[i] = 1.0 / std::max(polar::polarity[i], polmin);
    }
    device_array::copyin(n, polarity, polar::polarity);
    device_array::copyin(n, thole, polar::thole);
    device_array::copyin(n, pdamp, polar::pdamp);
    device_array::copyin(n, polarity_inv, pinvbuf.data());
  }
}

void dfield(real* gpu_field, real* gpu_fieldp) {
  if (epolar_electyp == elec_t::ewald)
    dfield_ewald(gpu_field, gpu_fieldp);
  else
    dfield_coulomb(gpu_field, gpu_fieldp);
}

void ufield(const real* gpu_uind, const real* gpu_uinp, real* gpu_field,
            real* gpu_fieldp) {
  if (epolar_electyp == elec_t::ewald)
    ufield_ewald(gpu_uind, gpu_uinp, gpu_field, gpu_fieldp);
  else
    ufield_coulomb(gpu_uind, gpu_uinp, gpu_field, gpu_fieldp);
}

void induce(real* gpu_ud, real* gpu_up) {
  induce_mutual_pcg1(gpu_ud, gpu_up);

  if (inform::debug && use_potent(polar_term)) {
    std::vector<double> uindbuf;
    uindbuf.resize(3 * n);
    DeviceMemory::copyout_array(uindbuf.data(), gpu_ud, 3 * n);
    bool header = true;
    for (int i = 0; i < n; ++i) {
      if (polar::polarity[i] != 0) {
        if (header) {
          header = false;
          print(stdout, "\n Induced Dipole Moments (Debye) :\n");
          print(stdout, "\n{0:4s}Atom{0:15s}X{0:12s}Y{0:12s}Z{0:11s}Total\n\n",
                "");
        }
        double u1 = uindbuf[3 * i];
        double u2 = uindbuf[3 * i + 1];
        double u3 = uindbuf[3 * i + 2];
        double unorm = std::sqrt(u1 * u1 + u2 * u2 + u3 * u3);
        u1 *= units::debye;
        u2 *= units::debye;
        u3 *= units::debye;
        unorm *= units::debye;
        print(stdout, "{:>8d}     {:13.4f}{:13.4f}{:13.4f} {:13.4f}\n", i + 1,
              u1, u2, u3, unorm);
      }
    }
  }
}

void epolar(int vers) {
  if (epolar_electyp == elec_t::coulomb)
    epolar_coulomb(vers);
  else if (epolar_electyp == elec_t::ewald)
    epolar_ewald(vers);
}
TINKER_NAMESPACE_END
