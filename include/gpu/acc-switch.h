#ifndef TINKER_GPU_ACC_SWITCH_H_
#define TINKER_GPU_ACC_SWITCH_H_

#include "acc-mathfunc.h"
#include "decl-real.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
enum {
  switch_default = 0,
  switch_vdw,
  switch_repuls,
  switch_disp,
  switch_charge,
  switch_chgdpl,
  switch_dipole,
  switch_mpole,
  switch_chgtrn,
  switch_ewald,
  switch_dewald,
  switch_usolve,
  switch_gkv,
  switch_gksa,
};

void switch_cut_off(int switch_type, double& cut, double& off);

#pragma acc routine seq
template <int DO_DTAPER>
void switch_taper5(real rik, real cut, real off, real& taper, real& dtaper) {
  // S2(x) = 6 x**5 - 15 x**4 + 10 x**3
  // S2(x): [0,1] :-> [0,1]
  // taper5: [cut,off] :-> [1,0]
  real _1_ab = REAL_RECIP(cut - off);
  real x = (rik - off) * _1_ab;
  real x2 = x * x;
  real x3 = x2 * x;
  taper = x3 * (6 * x2 - 15 * x + 10);
  if_constexpr(DO_DTAPER) { dtaper = 30 * REAL_SQ(x * (1 - x)) * _1_ab; }
}
}
TINKER_NAMESPACE_END

#endif
