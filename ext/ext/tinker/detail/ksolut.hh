#pragma once

#include "macro.h"

namespace tinker { namespace ksolut {
extern double*& solrad;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double* TINKER_MOD(ksolut, solrad);

double*& solrad = TINKER_MOD(ksolut, solrad);
#endif
} }
