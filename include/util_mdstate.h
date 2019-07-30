#ifndef TINKER_UTIL_MDSTATE_H_
#define TINKER_UTIL_MDSTATE_H_

#include "mod_mdstate.h"
#include "util_array.h"
#include "util_math.h"
#include "util_rc_man.h"

TINKER_NAMESPACE_BEGIN
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
TINKER_NAMESPACE_END

#include "mod_box.h"
#include "mod_couple.h"
#include "mod_nblist.h"
#include "mod_polgrp.h"
#include "util_random.h"

TINKER_NAMESPACE_BEGIN
void nblist_data(rc_t rc);
void goto_frame0(int idx0);
void goto_frame1(int idx1);
TINKER_NAMESPACE_END

#endif