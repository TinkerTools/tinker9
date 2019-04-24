#ifndef TINKER_MOD_KSTBND_HH_
#define TINKER_MOD_KSTBND_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kstbnd {
const int maxnsb = 2000;
extern double (&stbn)[maxnsb][2];
extern char (&ksb)[maxnsb][12];

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(kstbnd, stbn)[maxnsb][2];
extern "C" char m_tinker_mod(kstbnd, ksb)[maxnsb][12];

double (&stbn)[maxnsb][2] = m_tinker_mod(kstbnd, stbn);
char (&ksb)[maxnsb][12] = m_tinker_mod(kstbnd, ksb);
#endif

} TINKER_NAMESPACE_END

#endif
