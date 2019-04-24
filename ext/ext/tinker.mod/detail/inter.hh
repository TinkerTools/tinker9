#ifndef TINKER_MOD_INTER_HH_
#define TINKER_MOD_INTER_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace inter {
extern double& einter;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(inter, einter);

double& einter = m_tinker_mod(inter, einter);
#endif

} TINKER_NAMESPACE_END

#endif
