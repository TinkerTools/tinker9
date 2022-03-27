#include "ff/switch.h"
#include "math/inc.h"
#include <tinker/detail/limits.hh>
#include <tinker/detail/nonpol.hh>

namespace tinker {
real switchCut(Switch mode)
{
   real cut;
   using namespace limits;
   switch (mode) {
   case SWITCH_VDW:
      cut = vdwtaper;
      break;
   case SWITCH_REPULS:
      cut = reptaper;
      break;
   case SWITCH_DISP:
      cut = disptaper;
      break;
   case SWITCH_CHARGE:
      cut = chgtaper;
      break;
   case SWITCH_CHGDPL:
      cut = std::sqrt(chgtaper * dpltaper);
      break;
   case SWITCH_DIPOLE:
      cut = dpltaper;
      break;
   case SWITCH_MPOLE:
      cut = mpoletaper;
      break;
   case SWITCH_CHGTRN:
      cut = ctrntaper;
      break;
   case SWITCH_EWALD:
      cut = ewaldcut;
      break;
   case SWITCH_DEWALD:
      cut = dewaldcut;
      break;
   case SWITCH_USOLVE:
      cut = usolvcut;
      break;
   case SWITCH_GKV:
      cut = nonpol::spcut;
      break;
   case SWITCH_GKSA:
      cut = nonpol::stoff;
      break;
   default:
      cut = minOf(vdwtaper, reptaper, disptaper, chgtaper, dpltaper, mpoletaper, ctrntaper);
      break;
   }
   return cut;
}

real switchOff(Switch mode)
{
   real off;
   using namespace limits;
   switch (mode) {
   case SWITCH_VDW:
      off = vdwcut;
      break;
   case SWITCH_REPULS:
      off = repcut;
      break;
   case SWITCH_DISP:
      off = dispcut;
      break;
   case SWITCH_CHARGE:
      off = chgcut;
      break;
   case SWITCH_CHGDPL:
      off = std::sqrt(chgcut * dplcut);
      break;
   case SWITCH_DIPOLE:
      off = dplcut;
      break;
   case SWITCH_MPOLE:
      off = mpolecut;
      break;
   case SWITCH_CHGTRN:
      off = ctrncut;
      break;
   case SWITCH_EWALD:
      off = ewaldcut;
      break;
   case SWITCH_DEWALD:
      off = dewaldcut;
      break;
   case SWITCH_USOLVE:
      off = usolvcut;
      break;
   case SWITCH_GKV:
      off = nonpol::spoff;
      break;
   case SWITCH_GKSA:
      off = nonpol::stcut;
      break;
   default:
      off = minOf(vdwcut, repcut, dispcut, chgcut, dplcut, mpolecut, ctrncut);
      break;
   }
   return off;
}
}
