#ifndef TINKER_SWITCH_H_
#define TINKER_SWITCH_H_

#include "macro.h"

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

/**
 * @param[in] mode
 * potential type
 *
 * @param[out] cut
 * distance at which switching of the potential begins
 *
 * @param[out] off
 * distance at which the potential energy goes to zero
 */
void switch_cut_off(switch_t mode, double& cut, double& off);
TINKER_NAMESPACE_END

#endif
