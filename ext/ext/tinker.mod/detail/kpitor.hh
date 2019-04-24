#ifndef TINKER_MOD_KPITOR_HH_
#define TINKER_MOD_KPITOR_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kpitor {
const int maxnpt = 500;
extern double (&ptcon)[maxnpt];
extern char (&kpt)[maxnpt][8];

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(kpitor, ptcon)[maxnpt];
extern "C" char m_tinker_mod(kpitor, kpt)[maxnpt][8];

double (&ptcon)[maxnpt] = m_tinker_mod(kpitor, ptcon);
char (&kpt)[maxnpt][8] = m_tinker_mod(kpitor, kpt);
#endif

} TINKER_NAMESPACE_END

#endif
