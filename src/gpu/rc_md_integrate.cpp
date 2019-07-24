#include "gpu/decl_md.h"
#include "util/fort_str.h"
#include <ext/tinker/tinker_mod.h>

TINKER_NAMESPACE_BEGIN
namespace gpu {
static void (*integrator_)(int, real);

void integrate_data(rc_t rc) {
  if (rc & rc_copyin) {
    if (bath::isothermal) {
      fstr_view th = bath::thermostat;
      if (th == "BERENDSEN")
        thermostat = thermo_berendsen;
      else if (th == "BUSSI")
        thermostat = thermo_bussi;
      else if (th == "ANDERSEN")
        thermostat = thermo_andersen;
      else if (th == "NOSE-HOOVER")
        thermostat = thermo_nose_hoover_chain;
      else
        assert(false);
    } else {
      thermostat = thermo_null;
    }

    if (bath::isobaric) {
      fstr_view br = bath::barostat;
      if (br == "BERENDSEN")
        barostat = baro_berendsen;
      else if (br == "BUSSI")
        barostat = baro_bussi;
      else if (br == "NOSE-HOOVER")
        barostat = baro_nose_hoover_chain;
      else if (br == "MONTECARLO")
        barostat = baro_montecarlo;
      else
        assert(false);
    } else {
      barostat = baro_null;
    }

    fstr_view itg = mdstuf::integrate;
    integrator_ = nullptr;
    if (itg == "VERLET") {
      integrator_ = velocity_verlet;
    } else if (itg == "STOCHASTIC") {
    } else if (itg == "BAOAB") {
    } else if (itg == "BUSSI") {
    } else if (itg == "NOSE-HOOVER") {
    } else if (itg == "GHMC") {
    } else if (itg == "RIGIDBODY") {
    } else if (itg == "RESPA") {
    } else {
      // beeman
    }
  }
}

//======================================================================

extern void kinetic_acc_impl_(real& temp);
void kinetic(real& temp) { kinetic_acc_impl_(temp); }

extern void thermo_bussi_acc_impl_(real dt, real temp);
void temper(real dt, real& temp) {
  if (thermostat == thermo_null) {
    kinetic(temp);
    return;
  }

  if (thermostat == thermo_bussi)
    thermo_bussi_acc_impl_(dt, temp);
  else
    assert(false);

  kinetic(temp);
}

extern void mdrest_acc_impl_(int istep);
void mdrest(int istep) { mdrest_acc_impl_(istep); }

//======================================================================
// propagate

extern void propagate_xyz_acc_impl_(real dt);
void propagate_xyz(real dt) { propagate_xyz_acc_impl_(dt); }

extern void propagate_velocity_acc_impl_(real dt);
void propagate_velocity(real dt) { propagate_velocity_acc_impl_(dt); }

void propagate(int nsteps, real dt_ps, void (*itg)(int, real)) {
  if (itg == nullptr)
    itg = integrator_;

  for (int istep = 1; istep <= nsteps; ++istep) {
    itg(istep, dt_ps);

    // mdstat
    if (istep % inform::iwrite == 0)
      mdsave_async(istep, dt_ps);
    mdrest(istep);
  }
  mdsave_synchronize();
}

void md_data(rc_t rc) {
  if ((use_md & use_data) == 0)
    return;

  integrate_data(rc);
  mdsave_data(rc);
}
}
TINKER_NAMESPACE_END
