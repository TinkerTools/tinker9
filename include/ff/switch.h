#pragma once
#include "precision.h"

namespace tinker {
enum Switch
{
   SWITCH_DEFAULT,
   SWITCH_VDW,
   SWITCH_REPULS,
   SWITCH_DISP,
   SWITCH_CHARGE,
   SWITCH_CHGDPL,
   SWITCH_DIPOLE,
   SWITCH_MPOLE,
   SWITCH_CHGTRN,
   SWITCH_EWALD,
   SWITCH_DEWALD,
   SWITCH_USOLVE,
   SWITCH_GKV,
   SWITCH_GKSA,
};

/// \ingroup ff
/// \return Distance at which switching of the potential begins.
real switchCut(Switch mode);

/// \ingroup ff
/// \return Distance at which the potential energy goes to zero.
real switchOff(Switch mode);
}
