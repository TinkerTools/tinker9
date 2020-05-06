#pragma once

#include "macro.h"

namespace tinker { namespace kchrge {
extern double*& chg;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(kchrge, chg);

double*& chg = TINKER_MOD(kchrge, chg);
#endif
} }
