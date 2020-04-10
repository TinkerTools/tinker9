#pragma once

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace kvdwpr {
const int maxnvp = 500;
extern double (&radpr)[maxnvp];
extern double (&epspr)[maxnvp];
extern char (&kvpr)[maxnvp][8];

#ifdef TINKER_MOD_CPP_
extern "C" double TINKER_MOD(kvdwpr, radpr)[maxnvp];
extern "C" double TINKER_MOD(kvdwpr, epspr)[maxnvp];
extern "C" char TINKER_MOD(kvdwpr, kvpr)[maxnvp][8];

double (&radpr)[maxnvp] = TINKER_MOD(kvdwpr, radpr);
double (&epspr)[maxnvp] = TINKER_MOD(kvdwpr, epspr);
char (&kvpr)[maxnvp][8] = TINKER_MOD(kvdwpr, kvpr);
#endif
} TINKER_NAMESPACE_END
