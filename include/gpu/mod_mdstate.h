#ifndef TINKER_GPU_MOD_MDSTATE_H_
#define TINKER_GPU_MOD_MDSTATE_H_

#include "util_cxx.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
TINKER_EXTERN int use_data;
// number of frames
TINKER_EXTERN int trajn;
// number of atoms
TINKER_EXTERN int n;
// x y z coordinates
TINKER_EXTERN real *trajx, *trajy, *trajz;
TINKER_EXTERN real *x, *y, *z;
// velocities
TINKER_EXTERN real *vx, *vy, *vz;
// atomic mass
TINKER_EXTERN real *mass, *massinv;

// total potential energy
TINKER_EXTERN real* esum;
// total potential energy and kinetic energy
TINKER_EXTERN real epot, eksum, ekin[3][3];
// total gradients
TINKER_EXTERN real *gx, *gy, *gz;
// total virial
TINKER_EXTERN real* vir;
}
TINKER_NAMESPACE_END

#endif
