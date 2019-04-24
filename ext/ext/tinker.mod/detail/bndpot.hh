#ifndef TINKER_MOD_BNDPOT_HH_
#define TINKER_MOD_BNDPOT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace bndpot {
extern double& cbnd;
extern double& qbnd;
extern double& bndunit;
extern char (&bndtyp)[8];

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(bndpot, cbnd);
extern "C" double m_tinker_mod(bndpot, qbnd);
extern "C" double m_tinker_mod(bndpot, bndunit);
extern "C" char m_tinker_mod(bndpot, bndtyp)[8];

double& cbnd = m_tinker_mod(bndpot, cbnd);
double& qbnd = m_tinker_mod(bndpot, qbnd);
double& bndunit = m_tinker_mod(bndpot, bndunit);
char (&bndtyp)[8] = m_tinker_mod(bndpot, bndtyp);
#endif

} TINKER_NAMESPACE_END

#endif
