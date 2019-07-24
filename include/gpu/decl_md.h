#ifndef TINKER_GPU_DECL_MD_H_
#define TINKER_GPU_DECL_MD_H_

#include "decl_mdstate.h"
#include "mod_md.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
//======================================================================
// integrator

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
