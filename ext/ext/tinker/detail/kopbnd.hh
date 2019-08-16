#ifndef TINKER_MOD_KOPBND_HH_
#define TINKER_MOD_KOPBND_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace kopbnd {
const int maxnopb = 500;
extern double (&opbn)[maxnopb];
extern char (&kopb)[maxnopb][16];

#ifdef TINKER_MOD_CPP_
extern "C" double TINKER_MOD(kopbnd, opbn)[maxnopb];
extern "C" char TINKER_MOD(kopbnd, kopb)[maxnopb][16];

double (&opbn)[maxnopb] = TINKER_MOD(kopbnd, opbn);
char (&kopb)[maxnopb][16] = TINKER_MOD(kopbnd, kopb);
#endif
} TINKER_NAMESPACE_END

#endif
