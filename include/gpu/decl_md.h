#ifndef TINKER_GPU_DECL_MD_H_
#define TINKER_GPU_DECL_MD_H_

#include "decl_mdstate.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
//======================================================================
// integrator

typedef enum {
  thermo_berendsen,
  thermo_bussi,
  thermo_andersen,
  thermo_nose_hoover_chain,
  thermo_null
} thermostat_t;
extern thermostat_t thermostat;

typedef enum {
  baro_berendsen,
  baro_bussi,
  baro_nose_hoover_chain,
  baro_montecarlo,
  baro_null
} barostat_t;
extern barostat_t barostat;

void integrate_data(rc_t);

void kinetic(real& temp);
void temper(real dt, real& temp);
void mdrest(int istep);

void propagate_xyz(real dt);
void propagate_velocity(real dt);
void propagate(int nsteps, real dt_ps, void (*itg)(int, real) = nullptr);

//======================================================================
// mdsave

void mdsave_data(rc_t);

void mdsave_async(int istep, real dt);
void mdsave_synchronize();

//======================================================================

void md_data(rc_t rc);
}
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
namespace gpu {
void velocity_verlet(int istep, real dt_ps);
}
TINKER_NAMESPACE_END

#endif
