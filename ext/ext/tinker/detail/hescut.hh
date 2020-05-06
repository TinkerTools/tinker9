#pragma once

#include "macro.h"

namespace tinker { namespace hescut {
extern double& hesscut;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(hescut, hesscut);

double& hesscut = TINKER_MOD(hescut, hesscut);
#endif
} }
