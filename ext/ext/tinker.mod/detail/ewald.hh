#ifndef TINKER_MOD_EWALD_HH_
#define TINKER_MOD_EWALD_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace ewald {
extern double& aewald;
extern double& aeewald;
extern double& apewald;
extern double& adewald;
extern char (&boundary)[7];

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(ewald, aewald);
extern "C" double m_tinker_mod(ewald, aeewald);
extern "C" double m_tinker_mod(ewald, apewald);
extern "C" double m_tinker_mod(ewald, adewald);
extern "C" char m_tinker_mod(ewald, boundary)[7];

double& aewald = m_tinker_mod(ewald, aewald);
double& aeewald = m_tinker_mod(ewald, aeewald);
double& apewald = m_tinker_mod(ewald, apewald);
double& adewald = m_tinker_mod(ewald, adewald);
char (&boundary)[7] = m_tinker_mod(ewald, boundary);
#endif

} TINKER_NAMESPACE_END

#endif
