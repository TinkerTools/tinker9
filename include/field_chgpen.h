#pragma once
#include "macro.h"

namespace tinker {
// electrostatic field due to permanent multipoles
// clang-format off
void dfield_chgpen                  (real (*field)[3]);
void dfield_chgpen_nonewald         (real (*field)[3]);
void dfield_chgpen_ewald            (real (*field)[3]);
void dfield_chgpen_ewald_recip_self (real (*field)[3]);
void dfield_chgpen_ewald_real       (real (*field)[3]);

void dfield_chgpen_nonewald_acc         (real (*field)[3]);
void dfield_chgpen_ewald_real_acc       (real (*field)[3]);
void dfield_chgpen_nonewald_cu          (real (*field)[3]);
void dfield_chgpen_ewald_real_cu        (real (*field)[3]);
// clang-format on

// mutual electrostatic field due to induced dipole moments
// clang-format off
// -Tu operator
void ufield_chgpen                  (const real (*uind)[3], real (*field)[3]);
void ufield_chgpen_nonewald         (const real (*uind)[3], real (*field)[3]);
void ufield_chgpen_ewald            (const real (*uind)[3], real (*field)[3]);
void ufield_chgpen_ewald_recip_self (const real (*uind)[3], real (*field)[3]);
void ufield_chgpen_ewald_real       (const real (*uind)[3], real (*field)[3]);

void ufield_chgpen_nonewald_acc         (const real (*uind)[3], real (*field)[3]);
void ufield_chgpen_ewald_recip_self_acc (const real (*uind)[3], real (*field)[3]);
void ufield_chgpen_ewald_real_acc       (const real (*uind)[3], real (*field)[3]);
void ufield_chgpen_nonewald_cu          (const real (*uind)[3], real (*field)[3]);
void ufield_chgpen_ewald_real_cu        (const real (*uind)[3], real (*field)[3]);
// clang-format on
}
