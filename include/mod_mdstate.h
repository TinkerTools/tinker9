#ifndef TINKER_MOD_MDSTATE_H_
#define TINKER_MOD_MDSTATE_H_

#include "util_cxx.h"

TINKER_NAMESPACE_BEGIN
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

constexpr int _x = 0; /// x direction
constexpr int _y = 1; /// y direction
constexpr int _z = 2; /// z direction

constexpr int _xx = 0; /// xx component
constexpr int _yx = 1; /// yx component
constexpr int _zx = 2; /// zx component
constexpr int _xy = 3; /// xy component
constexpr int _yy = 4; /// yy component
constexpr int _zy = 5; /// zy component
constexpr int _xz = 6; /// xz component
constexpr int _yz = 7; /// yz component
constexpr int _zz = 8; /// zz component
TINKER_NAMESPACE_END

#endif
