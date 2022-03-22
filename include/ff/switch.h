#pragma once
#include "macro.h"

namespace tinker {
enum switch_t
{
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
};

/// \return Distance at which switching of the potential begins.
real switch_cut(switch_t mode);

/// \return Distance at which the potential energy goes to zero.
real switch_off(switch_t mode);
}
