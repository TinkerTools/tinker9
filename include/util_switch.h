#ifndef TINKER_UTIL_SWITCH_H_
#define TINKER_UTIL_SWITCH_H_

#include "util_macro.h"

TINKER_NAMESPACE_BEGIN
typedef enum {
  switch_default,
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
} switch_t;

void switch_cut_off(switch_t switch_type, double& cut, double& off);
TINKER_NAMESPACE_END

#endif
