#pragma once

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace charge {
extern int& nion;
extern int*& iion;
extern int*& jion;
extern int*& kion;
extern double*& pchg;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(charge, nion);
extern "C" int* TINKER_MOD(charge, iion);
extern "C" int* TINKER_MOD(charge, jion);
extern "C" int* TINKER_MOD(charge, kion);
extern "C" double* TINKER_MOD(charge, pchg);

int& nion = TINKER_MOD(charge, nion);
int*& iion = TINKER_MOD(charge, iion);
int*& jion = TINKER_MOD(charge, jion);
int*& kion = TINKER_MOD(charge, kion);
double*& pchg = TINKER_MOD(charge, pchg);
#endif
} TINKER_NAMESPACE_END
