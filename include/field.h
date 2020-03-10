#pragma once
#include "macro.h"


TINKER_NAMESPACE_BEGIN
// electrostatic field due to permanent multipoles
// clang-format off
void dfield                  (real (*field)[3], real (*fieldp)[3]);
void dfield_nonewald         (real (*field)[3], real (*fieldp)[3]);
void dfield_ewald            (real (*field)[3], real (*fieldp)[3]);
void dfield_ewald_recip_self (real (*field)[3], real (*fieldp)[3]);
void dfield_ewald_real       (real (*field)[3], real (*fieldp)[3]);

void dfield_nonewald_acc         (real (*field)[3], real (*fieldp)[3]);
void dfield_ewald_recip_self_acc (real (*field)[3]);
void dfield_ewald_real_acc       (real (*field)[3], real (*fieldp)[3]);
void dfield_nonewald_cu          (real (*field)[3], real (*fieldp)[3]);
void dfield_ewald_real_cu        (real (*field)[3], real (*fieldp)[3]);
// clang-format on


// mutual electrostatic field due to induced dipole moments
// clang-format off
// -Tu operator
void ufield                  (const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufield_nonewald         (const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufield_ewald            (const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufield_ewald_recip_self (const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufield_ewald_real       (const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);

void ufield_nonewald_acc         (const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufield_ewald_recip_self_acc (const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufield_ewald_real_acc       (const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufield_nonewald_cu          (const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufield_ewald_real_cu        (const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
// clang-format on
TINKER_NAMESPACE_END
