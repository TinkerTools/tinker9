#ifndef TINKER_MOD_KPITOR_HH_
#define TINKER_MOD_KPITOR_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace kpitor {
const int maxnpt = 500;
extern double (&ptcon)[maxnpt];
extern char (&kpt)[maxnpt][8];

#ifdef TINKER_MOD_CPP_
extern "C" double TINKER_MOD(kpitor, ptcon)[maxnpt];
extern "C" char TINKER_MOD(kpitor, kpt)[maxnpt][8];

double (&ptcon)[maxnpt] = TINKER_MOD(kpitor, ptcon);
char (&kpt)[maxnpt][8] = TINKER_MOD(kpitor, kpt);
#endif
} TINKER_NAMESPACE_END

#endif
