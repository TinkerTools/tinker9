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

TINKER_EXTERN thermostat_t thermostat;
TINKER_EXTERN barostat_t barostat;

enum {
  _x = 0, /// x direction
  _y = 1, /// y direction
  _z = 2, /// z direction

  _xx = 0, /// xx component
  _yx = 1, /// yx component
  _zx = 2, /// zx component
  _xy = 3, /// xy component
  _yy = 4, /// yy component
  _zy = 5, /// zy component
  _xz = 6, /// xz component
  _yz = 7, /// yz component
  _zz = 8  /// zz component
};

namespace calc {
enum {
  xyz = 0x001,  /// xyz
  vel = 0x002,  /// velocity
  mass = 0x004, /// mass
  traj = 0x008, /// trajectory

  energy = 0x010, /// energy 16
  grad = 0x020,   /// gradient 32
  virial = 0x040, /// virial 64
  analyz = 0x080, /// analyze 128

  md = 0x100,

  // clang-format off
  vmask = energy + grad + virial +  analyz,
  v0 = energy,                 ///  16
  v1 = energy + grad + virial, /// 112
  v3 = energy + analyz,        /// 144
  v4 = energy + grad,          ///  48
  v5 = grad,                   ///  32
  v6 = grad + virial,          ///  96
  // clang-format on
};
}
TINKER_NAMESPACE_END

#endif
