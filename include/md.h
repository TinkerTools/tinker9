#ifndef TINKER_MD_H_
#define TINKER_MD_H_

#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
TINKER_EXTERN int rc_flag;
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
  xyz = 0x001,  ///< xyz
  vel = 0x002,  ///< velocity
  mass = 0x004, ///< mass
  traj = 0x008, ///< trajectory

  energy = 0x010, ///< energy 16
  grad = 0x020,   ///< gradient 32
  virial = 0x040, ///< virial 64
  analyz = 0x080, ///< analyze 128

  md = 0x100, ///< md calculation

  vmask = energy + grad + virial + analyz,
  v0 = energy,                 //  16
  v1 = energy + grad + virial, // 112
  v3 = energy + analyz,        // 144
  v4 = energy + grad,          //  48
  v5 = grad,                   //  32
  v6 = grad + virial,          //  96
};
}

template <int USE>
void sanity_check() {
  constexpr int do_e = USE & calc::energy;
  constexpr int do_a = USE & calc::analyz;
  constexpr int do_g = USE & calc::grad;
  constexpr int do_v = USE & calc::virial;
  // if calc::virial, must calc::grad
  static_assert(do_v ? do_g : true, "");
  // if calc::analyz, must calc::energy
  static_assert(do_a ? do_e : true, "");
}

//======================================================================
void egv_data(rc_op op);

//======================================================================
// energy, gradient, and virial de/allocation

void alloc_ev(real** gpu_e, real** gpu_v);
void free_ev(real* gpu_e, real* gpu_v);

void alloc_nev(int** gpu_ne, real** gpu_e, real** gpu_v);
void free_nev(int* gpu_ne, real* gpu_e, real* gpu_v);

double get_energy(const real* e_gpu);
int get_count(const int* ecount_gpu);
void get_virial(double* v_out, const real* v_gpu);
/// zero out global total energy, gradients, and virial on device
void zero_egv(int vers);
void zero_egv();

/// sum potential energies and virials
void sum_energies(int vers);

void goto_frame0(int idx0);
void goto_frame1(int idx1);
TINKER_NAMESPACE_END

#endif
