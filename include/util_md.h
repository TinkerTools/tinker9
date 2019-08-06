#ifndef TINKER_UTIL_MD_H_
#define TINKER_UTIL_MD_H_

#include "rc_man.h"
#include <string>

TINKER_NAMESPACE_BEGIN
//======================================================================
// integrator

void integrate_data(rc_op);

void kinetic(real& temp);
void temper(real dt, real& temp);
void mdrest(int istep);

void propagate_xyz(real dt);
void propagate_velocity(real dt);
void propagate(int nsteps, real dt_ps, void (*itg)(int, real) = nullptr);

//======================================================================
// mdsave

void mdsave_data(rc_op);

void mdsave_async(int istep, real dt);
void mdsave_synchronize();

//======================================================================

void md_data(rc_op op);

void goto_frame(int idx0);
void copyin_tinker_arc(const std::string& arcfile, int first1, int last1,
                       int step);
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
void velocity_verlet(int istep, real dt_ps);
TINKER_NAMESPACE_END

#endif
