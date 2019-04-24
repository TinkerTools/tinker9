#ifndef TINKER_MOD_CHRONO_HH_
#define TINKER_MOD_CHRONO_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace chrono {
extern double& twall;
extern double& tcpu;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(chrono, twall);
extern "C" double m_tinker_mod(chrono, tcpu);

double& twall = m_tinker_mod(chrono, twall);
double& tcpu = m_tinker_mod(chrono, tcpu);
#endif

} TINKER_NAMESPACE_END

#endif
