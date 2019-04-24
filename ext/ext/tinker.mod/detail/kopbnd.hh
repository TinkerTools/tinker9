#ifndef TINKER_MOD_KOPBND_HH_
#define TINKER_MOD_KOPBND_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kopbnd {
const int maxnopb = 500;
extern double (&opbn)[maxnopb];
extern char (&kopb)[maxnopb][16];

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(kopbnd, opbn)[maxnopb];
extern "C" char m_tinker_mod(kopbnd, kopb)[maxnopb][16];

double (&opbn)[maxnopb] = m_tinker_mod(kopbnd, opbn);
char (&kopb)[maxnopb][16] = m_tinker_mod(kopbnd, kopb);
#endif

} TINKER_NAMESPACE_END

#endif
