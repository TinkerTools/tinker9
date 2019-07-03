#ifndef TINKER_GPU_DECL_MDSTATE_H_
#define TINKER_GPU_DECL_MDSTATE_H_

#include "decl_basic.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
enum {
  use_xyz = 0x001,   /// xyz
  use_vel = 0x002,   /// velocity
  use_accel = 0x004, /// acceleration
  use_mass = 0x008,  /// mass

  use_energy = 0x010, /// energy 16
  use_grad = 0x020,   /// gradient 32
  use_virial = 0x040, /// virial 64
  use_analyz = 0x080, /// analyze 128

  // clang-format off
  vmask = use_energy + use_grad + use_virial + use_analyz,
  v0 = use_energy,                         ///  16
  v1 = use_energy + use_grad + use_virial, /// 112
  v3 = use_energy + use_analyz,            /// 144
  v4 = use_energy + use_grad,              ///  48
  v5 = use_grad,                           ///  32
  v6 = use_grad + use_virial,              ///  96
  // clang-format on

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

template <int USE>
void sanity_check() {
  constexpr int do_e = USE & use_energy;
  constexpr int do_a = USE & use_analyz;
  constexpr int do_g = USE & use_grad;
  constexpr int do_v = USE & use_virial;
  // if use_virial, must use_grad
  static_assert(do_v ? do_g : true, "");
  // if use_analyz, must use_energy
  static_assert(do_a ? do_e : true, "");
}

extern int use_data;

//======================================================================
// number of atoms

extern int n;
void n_data(rc_t rc);

//======================================================================
/// x y z coordinates
extern real *x, *y, *z;
void xyz_data(rc_t rc);

//======================================================================
/// velocities
extern real *vx, *vy, *vz;
void vel_data(rc_t rc);

//======================================================================
/// accelerations
extern real *ax, *ay, *az;
void accel_data(rc_t rc);

//======================================================================
/// atomic mass
extern real* mass;
extern real* massinv;
void mass_data(rc_t rc);

//======================================================================
/// total potential energy
extern real* esum;
/// total gradients
extern real *gx, *gy, *gz;
/// total virial
extern real* vir;

void egv_data(rc_t rc);

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
}
TINKER_NAMESPACE_END

#include "decl_box.h"
#include "decl_couple.h"
#include "decl_nblist.h"
#include "decl_polgrp.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void mdstate_data(rc_t rc);
}
TINKER_NAMESPACE_END

#endif
