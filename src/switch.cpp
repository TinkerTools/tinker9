#include "switch.h"
#include "mathfunc.h"
#include <tinker/detail/limits.hh>
#include <tinker/detail/nonpol.hh>


TINKER_NAMESPACE_BEGIN
real switch_cut(switch_t mode)
{
   real cut;
   using namespace limits;
   switch (mode) {
   case switch_vdw:
      cut = vdwtaper;
      break;
   case switch_repuls:
      cut = reptaper;
      break;
   case switch_disp:
      cut = disptaper;
      break;
   case switch_charge:
      cut = chgtaper;
      break;
   case switch_chgdpl:
      cut = std::sqrt(chgtaper * dpltaper);
      break;
   case switch_dipole:
      cut = dpltaper;
      break;
   case switch_mpole:
      cut = mpoletaper;
      break;
   case switch_chgtrn:
      cut = ctrntaper;
      break;
   case switch_ewald:
      cut = ewaldcut;
      break;
   case switch_dewald:
      cut = dewaldcut;
      break;
   case switch_usolve:
      cut = usolvcut;
      break;
   case switch_gkv:
      cut = nonpol::spcut;
      break;
   case switch_gksa:
      cut = nonpol::stoff;
      break;
   default: // switch_default
      cut = min_of(vdwtaper, reptaper, disptaper, chgtaper, dpltaper,
                   mpoletaper, ctrntaper);
      break;
   }
   return cut;
}


real switch_off(switch_t mode)
{
   real off;
   using namespace limits;
   switch (mode) {
   case switch_vdw:
      off = vdwcut;
      break;
   case switch_repuls:
      off = repcut;
      break;
   case switch_disp:
      off = dispcut;
      break;
   case switch_charge:
      off = chgcut;
      break;
   case switch_chgdpl:
      off = std::sqrt(chgcut * dplcut);
      break;
   case switch_dipole:
      off = dplcut;
      break;
   case switch_mpole:
      off = mpolecut;
      break;
   case switch_chgtrn:
      off = ctrncut;
      break;
   case switch_ewald:
      off = ewaldcut;
      break;
   case switch_dewald:
      off = dewaldcut;
      break;
   case switch_usolve:
      off = usolvcut;
      break;
   case switch_gkv:
      off = nonpol::spoff;
      break;
   case switch_gksa:
      off = nonpol::stcut;
      break;
   default: // switch_default
      off = min_of(vdwcut, repcut, dispcut, chgcut, dplcut, mpolecut, ctrncut);
      break;
   }
   return off;
}
TINKER_NAMESPACE_END
