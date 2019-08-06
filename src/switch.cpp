#include "mathfunc.h"
#include "switch.h"
#include <ext/tinker/tinker_mod.h>

TINKER_NAMESPACE_BEGIN
void switch_cut_off(switch_t switch_type, double& cut, double& off) {
  using namespace limits;
  switch (switch_type) {
  case switch_vdw:
    off = vdwcut;
    cut = vdwtaper;
    break;
  case switch_repuls:
    off = repcut;
    cut = reptaper;
    break;
  case switch_disp:
    off = dispcut;
    cut = disptaper;
    break;
  case switch_charge:
    off = chgcut;
    cut = chgtaper;
    break;
  case switch_chgdpl:
    off = std::sqrt(chgcut * dplcut);
    cut = std::sqrt(chgtaper * dpltaper);
    break;
  case switch_dipole:
    off = dplcut;
    cut = dpltaper;
    break;
  case switch_mpole:
    off = mpolecut;
    cut = mpoletaper;
    break;
  case switch_chgtrn:
    off = ctrncut;
    cut = ctrntaper;
    break;
  case switch_ewald:
    off = ewaldcut;
    cut = ewaldcut;
    break;
  case switch_dewald:
    off = dewaldcut;
    cut = dewaldcut;
    break;
  case switch_usolve:
    off = usolvcut;
    cut = usolvcut;
    break;
  case switch_gkv:
    off = nonpol::spoff;
    cut = nonpol::spcut;
    break;
  case switch_gksa:
    off = nonpol::stcut;
    cut = nonpol::stoff;
    break;
  default: // switch_default
    off = min_of(vdwcut, repcut, dispcut, chgcut, dplcut, mpolecut, ctrncut);
    cut = min_of(vdwtaper, reptaper, disptaper, chgtaper, dpltaper, mpoletaper,
                 ctrntaper);
    break;
  }
}
TINKER_NAMESPACE_END
