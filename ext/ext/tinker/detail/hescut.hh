#pragma once

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace hescut {
extern double& hesscut;

#ifdef TINKER_MOD_CPP_
extern "C" double TINKER_MOD(hescut, hesscut);

double& hesscut = TINKER_MOD(hescut, hesscut);
#endif
} TINKER_NAMESPACE_END
