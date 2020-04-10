#pragma once

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace chrono {
extern double& twall;
extern double& tcpu;

#ifdef TINKER_MOD_CPP_
extern "C" double TINKER_MOD(chrono, twall);
extern "C" double TINKER_MOD(chrono, tcpu);

double& twall = TINKER_MOD(chrono, twall);
double& tcpu = TINKER_MOD(chrono, tcpu);
#endif
} TINKER_NAMESPACE_END
