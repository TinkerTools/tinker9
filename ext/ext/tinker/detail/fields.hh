#ifndef TINKER_MOD_FIELDS_HH_
#define TINKER_MOD_FIELDS_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace fields {
const int maxbio = 10000;
extern int (&biotyp)[maxbio];
extern char (&forcefield)[20];

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(fields, biotyp)[maxbio];
extern "C" char TINKER_MOD(fields, forcefield)[20];

int (&biotyp)[maxbio] = TINKER_MOD(fields, biotyp);
char (&forcefield)[20] = TINKER_MOD(fields, forcefield);
#endif
} TINKER_NAMESPACE_END

#endif
