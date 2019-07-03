#include "gpu/decl_md.h"
#include "gpu/decl_mdstate.h"
#include "gpu/decl_potent.h"
#include "gpu/e_potential.h"
#include "gpu/rc.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void mdsave(int istep, real dt, real epot, real eksum) {
  int modsave = istep % inform::iwrite;
  if (modsave != 0)
    return;

  if (bound::use_bounds)
    box_data(rc_copyout);

  // TODO: RIGIDBODY
  // TODO: rc_copyout_async
  const rc_t rc = rc_copyout;

  xyz_data(rc);
  if (output::velsave)
    vel_data(rc);

  // if (output::frcsave)
  // TODO: frcsave

  if (output::uindsave && use_potent(polar_term))
    copyout_array(polar::uind, &uind[0][0], 3 * n);

  double dt1 = dt;
  double epot1 = epot;
  double eksum1 = eksum;
  TINKER_RT(mdsave)(&istep, &dt1, &epot1, &eksum1);
}

extern void mdrest_acc_impl__(int istep);
void mdrest(int istep) { mdrest_acc_impl__(istep); }

extern void propagate_xyz_acc_impl__(real dt);
void propagate_xyz(real dt) { propagate_xyz_acc_impl__(dt); }

extern void propagate_velocity_acc_impl__(real dt);
void propagate_velocity(real dt) { propagate_velocity_acc_impl__(dt); }

// integrators

void velocity_verlet(int istep, real dt) {
  real dt_2 = 0.5f * dt;

  // v += a * dt/2
  propagate_velocity(dt_2);

  // s += v * dt
  propagate_xyz(dt);
  // update gradient
  energy_potential(use_data & vmask);

  // v += a * dt/2
  propagate_velocity(dt_2);

  // mdstat

  // mdsave

  mdrest(istep);
}
}
TINKER_NAMESPACE_END
