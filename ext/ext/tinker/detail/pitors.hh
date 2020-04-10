#pragma once

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace pitors {
extern int& npitors;
extern int*& ipit;
extern double*& kpit;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(pitors, npitors);
extern "C" int* TINKER_MOD(pitors, ipit);
extern "C" double* TINKER_MOD(pitors, kpit);

int& npitors = TINKER_MOD(pitors, npitors);
int*& ipit = TINKER_MOD(pitors, ipit);
double*& kpit = TINKER_MOD(pitors, kpit);
#endif
} TINKER_NAMESPACE_END
