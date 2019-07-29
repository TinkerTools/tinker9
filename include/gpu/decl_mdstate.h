#ifndef TINKER_GPU_DECL_MDSTATE_H_
#define TINKER_GPU_DECL_MDSTATE_H_

#include "gpu_array.h"
#include "mod_mdstate.h"
#include "util_math.h"
#include "util_rc_man.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
enum {
  use_xyz = 0x001,  /// xyz
  use_vel = 0x002,  /// velocity
  use_mass = 0x004, /// mass
  use_traj = 0x008, /// trajectory

  use_energy = 0x010, /// energy 16
  use_grad = 0x020,   /// gradient 32
  use_virial = 0x040, /// virial 64
  use_analyz = 0x080, /// analyze 128

  use_md = 0x100,

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

//======================================================================
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

#include "mod_box.h"
#include "mod_couple.h"
#include "mod_nblist.h"
#include "mod_polgrp.h"
#include "util_random.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void nblist_data(rc_t rc);
void goto_frame0(int idx0);
void goto_frame1(int idx1);
}
TINKER_NAMESPACE_END

#endif
