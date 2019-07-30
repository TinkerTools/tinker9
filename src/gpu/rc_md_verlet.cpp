#include "gpu/e_potential.h"
#include "util_md.h"
#include "util_mdstate.h"
#include <ext/tinker/tinker_mod.h>

TINKER_NAMESPACE_BEGIN
void velocity_verlet(int istep, real dt_ps) {
  bool save = !(istep % inform::iwrite);
  int vers0 = use_data & (calc::virial | calc::grad | calc::energy);
  int vers1 = use_data & (calc::virial | calc::grad);

  real dt_2 = 0.5f * dt_ps;

  // v += a * dt/2
  propagate_velocity(dt_2);

  // s += v * dt
  propagate_xyz(dt_ps);
  // update gradient
  if (save)
    energy_potential(vers0);
  else
    energy_potential(vers1);

  // v += a * dt/2
  propagate_velocity(dt_2);

  real temp;
  temper(dt_ps, temp);
}
TINKER_NAMESPACE_END
