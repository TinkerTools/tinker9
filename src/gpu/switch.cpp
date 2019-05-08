#include "gpu/switch.h"
#include "tinker.mod.h"
#include "util/math.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
// cut: taper
// off: cutoff
// c: at least c[6]
static double switch_c05_no_check_(double cut, double off, double* c) {
  c[0] = 0;
  c[1] = 0;
  c[2] = 0;
  c[3] = 0;
  c[4] = 0;
  c[5] = 0;
  // if (cut < off) {
  double off2 = REAL_SQ(off);
  double cut2 = REAL_SQ(cut);
  double _1_denom = REAL_POW(off - cut, -5);

  c[0] = off * off2 * (off2 - 5 * off * cut + 10 * cut2) * _1_denom;
  c[1] = -30 * off2 * cut2 * _1_denom;
  c[2] = 30 * (off2 * cut + off * cut2) * _1_denom;
  c[3] = -10 * (off2 + 4 * off * cut + cut2) * _1_denom;
  c[4] = 15 * (off + cut) * _1_denom;
  c[5] = -6 * _1_denom;
  // }
  return cut2;
}

// cut: taper
// off: cutoff
// f: at least f[8]
static double switch_f07_no_check_(double cut, double off, double* f) {
  f[0] = 0;
  f[1] = 0;
  f[2] = 0;
  f[3] = 0;
  f[4] = 0;
  f[5] = 0;
  f[6] = 0;
  f[7] = 0;
  // if (cut < off) {
  double off2 = off * off;
  double off3 = off2 * off;
  double off4 = off2 * off2;
  double off5 = off2 * off3;
  double off6 = off3 * off3;
  double off7 = off3 * off4;
  double cut2 = cut * cut;
  double cut3 = cut2 * cut;
  double cut4 = cut2 * cut2;
  double cut5 = cut2 * cut3;
  double cut6 = cut3 * cut3;
  double cut7 = cut3 * cut4;
  double term = 9.3 * cut * off / (off - cut);
  double denom = cut7 - 7 * cut6 * off + 21 * cut5 * off2 - 35 * cut4 * off3 +
      35 * cut3 * off4 - 21 * cut2 * off5 + 7 * cut * off6 - off7;
  denom *= term;
  denom = REAL_RECIP(denom);

  f[0] = cut3 * off3 * (-39 * cut + 64 * off) * denom;
  f[1] = cut2 * off2 * (117 * cut2 - 100 * cut * off - 192 * off2) * denom;
  f[2] = cut * off *
      (-117 * cut3 - 84 * cut2 * off + 534 * cut * off2 + 192 * off3) * denom;
  f[3] = (39 * cut4 + 212 * cut3 * off - 450 * cut2 * off2 - 612 * cut * off3 -
          64 * off4) *
      denom;
  f[4] = (-92 * cut3 + 66 * cut2 * off + 684 * cut * off2 + 217 * off3) * denom;
  f[5] = (42 * cut2 - 300 * cut * off - 267 * off2) * denom;
  f[6] = (36 * cut + 139 * off) * denom;
  f[7] = -25 * denom;
  // }
  return cut2;
}

void switching(int switch_type, double* coeff, double& taper2) {
  using namespace limits;
  double off, cut;

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
  default /* switch_default */:
    off = min_of(vdwcut, repcut, dispcut, chgcut, dplcut, mpolecut, ctrncut);
    cut = min_of(vdwtaper, reptaper, disptaper, chgtaper, dpltaper, mpoletaper,
                 ctrntaper);
    break;
  }

  if (cut < off) {
    if (switch_type != switch_charge)
      taper2 = switch_c05_no_check_(cut, off, coeff);
    else
      taper2 = switch_f07_no_check_(cut, off, coeff);
  }
}
}
TINKER_NAMESPACE_END
