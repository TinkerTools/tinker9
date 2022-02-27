#pragma once
#include "macro.h"


namespace tinker {
// electrostatic field due to permanent multipoles
// clang-format off
void dfield_aplus                  (real (*field)[3]);
void dfield_aplus_nonewald         (real (*field)[3]);
void dfield_aplus_ewald            (real (*field)[3]);
void dfield_aplus_ewald_recip_self (real (*field)[3]);
void dfield_aplus_ewald_real       (real (*field)[3]);

void dfield_aplus_nonewald_acc         (real (*field)[3]);
void dfield_aplus_ewald_real_acc       (real (*field)[3]);
void dfield_aplus_nonewald_cu          (real (*field)[3]);
void dfield_aplus_ewald_real_cu        (real (*field)[3]);
// clang-format on


// mutual electrostatic field due to induced dipole moments
// clang-format off
// -Tu operator
void ufield_aplus                  (const real (*uind)[3], real (*field)[3]);
void ufield_aplus_nonewald         (const real (*uind)[3], real (*field)[3]);
void ufield_aplus_ewald            (const real (*uind)[3], real (*field)[3]);
void ufield_aplus_ewald_recip_self (const real (*uind)[3], real (*field)[3]);
void ufield_aplus_ewald_real       (const real (*uind)[3], real (*field)[3]);

void ufield_aplus_nonewald_acc         (const real (*uind)[3], real (*field)[3]);
void ufield_aplus_ewald_recip_self_acc (const real (*uind)[3], real (*field)[3]);
void ufield_aplus_ewald_real_acc       (const real (*uind)[3], real (*field)[3]);
void ufield_aplus_nonewald_cu          (const real (*uind)[3], real (*field)[3]);
void ufield_aplus_ewald_real_cu        (const real (*uind)[3], real (*field)[3]);
// clang-format on
}
