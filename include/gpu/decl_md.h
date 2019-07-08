#ifndef TINKER_GPU_DECL_MD_H_
#define TINKER_GPU_DECL_MD_H_

#include "decl_basic.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
typedef enum {
  thermo_berendsen,
  thermo_bussi,
  thermo_andersen,
  thermo_nose_hoover_chain,
  thermo_null
} thermostat_t;

typedef enum {
  baro_berendsen,
  baro_bussi,
  baro_nose_hoover_chain,
  baro_montecarlo,
  baro_null
} barostat_t;

void mdinit();
/// pass device data back to Tinker mdsave subroutine
void mdsave(int istep, real dt, real epot, real eksum);

/**
 * finds and removes any translational or rotational kinetic energy of the
 * overall system center of mass
 */
void mdrest(int istep);

void propagate_xyz(real dt);
void propagate_velocity(real dt);

// integrators

void velocity_verlet(int istep, real dt);
}
TINKER_NAMESPACE_END

#endif
