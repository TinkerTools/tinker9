#ifndef TINKER_MOD_URYPOT_HH_
#define TINKER_MOD_URYPOT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace urypot {
extern double& cury;
extern double& qury;
extern double& ureyunit;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(urypot, cury);
extern "C" double m_tinker_mod(urypot, qury);
extern "C" double m_tinker_mod(urypot, ureyunit);

double& cury = m_tinker_mod(urypot, cury);
double& qury = m_tinker_mod(urypot, qury);
double& ureyunit = m_tinker_mod(urypot, ureyunit);
#endif

} TINKER_NAMESPACE_END

#endif
