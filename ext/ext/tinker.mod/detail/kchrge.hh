#ifndef TINKER_MOD_KCHRGE_HH_
#define TINKER_MOD_KCHRGE_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kchrge {
extern double*& chg;

#ifdef TINKER_MOD_CPP_
extern "C" double* m_tinker_mod(kchrge, chg);

double*& chg = m_tinker_mod(kchrge, chg);
#endif

} TINKER_NAMESPACE_END

#endif
