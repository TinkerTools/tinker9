#pragma once

#include "macro.h"

namespace tinker { namespace fields {
const int maxbio = 10000;
extern int (&biotyp)[maxbio];
extern char (&forcefield)[20];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int TINKER_MOD(fields, biotyp)[maxbio];
extern "C" char TINKER_MOD(fields, forcefield)[20];

int (&biotyp)[maxbio] = TINKER_MOD(fields, biotyp);
char (&forcefield)[20] = TINKER_MOD(fields, forcefield);
#endif
} }
