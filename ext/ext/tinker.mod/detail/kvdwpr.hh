#ifndef TINKER_MOD_KVDWPR_HH_
#define TINKER_MOD_KVDWPR_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kvdwpr {
const int maxnvp = 500;
extern double (&radpr)[maxnvp];
extern double (&epspr)[maxnvp];
extern char (&kvpr)[maxnvp][8];

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(kvdwpr, radpr)[maxnvp];
extern "C" double m_tinker_mod(kvdwpr, epspr)[maxnvp];
extern "C" char m_tinker_mod(kvdwpr, kvpr)[maxnvp][8];

double (&radpr)[maxnvp] = m_tinker_mod(kvdwpr, radpr);
double (&epspr)[maxnvp] = m_tinker_mod(kvdwpr, epspr);
char (&kvpr)[maxnvp][8] = m_tinker_mod(kvdwpr, kvpr);
#endif

} TINKER_NAMESPACE_END

#endif
