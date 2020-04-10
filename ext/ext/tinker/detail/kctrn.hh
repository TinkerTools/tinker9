#pragma once

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace kctrn {
extern double*& ctchg;
extern double*& ctdmp;

#ifdef TINKER_MOD_CPP_
extern "C" double* TINKER_MOD(kctrn, ctchg);
extern "C" double* TINKER_MOD(kctrn, ctdmp);

double*& ctchg = TINKER_MOD(kctrn, ctchg);
double*& ctdmp = TINKER_MOD(kctrn, ctdmp);
#endif
} TINKER_NAMESPACE_END
