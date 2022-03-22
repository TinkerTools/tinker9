#pragma once
#include "ff/molecule.h"

namespace tinker {
const int couple_maxn12 = 8;
TINKER_EXTERN int (*couple_i12)[couple_maxn12];
TINKER_EXTERN int* couple_n12;

TINKER_EXTERN Molecule molecule;

TINKER_EXTERN Group grp;
}
