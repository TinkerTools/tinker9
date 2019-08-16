#ifndef TINKER_MOD_KCHRGE_HH_
#define TINKER_MOD_KCHRGE_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace kchrge {
extern double*& chg;

#ifdef TINKER_MOD_CPP_
extern "C" double* TINKER_MOD(kchrge, chg);

double*& chg = TINKER_MOD(kchrge, chg);
#endif
} TINKER_NAMESPACE_END

#endif
