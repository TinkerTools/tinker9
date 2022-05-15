#pragma once
#include "ff/precision.h"

namespace tinker {
/// \ingroup hippopolar
/// \{
// electrostatic field due to permanent multipoles
void dfieldChgpen(real (*field)[3]);

// mutual electrostatic field due to induced dipole moments
// -Tu operator
void ufieldChgpen(const real (*uind)[3], real (*field)[3]);
void ufieldChgpenEwaldRecipSelf(const real (*uind)[3], real (*field)[3]);

void diagPrecond2(const real (*rsd)[3], real (*zrsd)[3]);

void sparsePrecondBuild2();
void sparsePrecondApply2(const real (*rsd)[3], real (*zrsd)[3]);

void induce2(real (*uind)[3]);
/// \}
}
