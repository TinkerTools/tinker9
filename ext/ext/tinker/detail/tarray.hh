#ifndef TINKER_MOD_TARRAY_HH_
#define TINKER_MOD_TARRAY_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace tarray {
extern int& ntpair;
extern int*& tindex;
extern double*& tdipdip;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(tarray, ntpair);
extern "C" int* TINKER_MOD(tarray, tindex);
extern "C" double* TINKER_MOD(tarray, tdipdip);

int& ntpair = TINKER_MOD(tarray, ntpair);
int*& tindex = TINKER_MOD(tarray, tindex);
double*& tdipdip = TINKER_MOD(tarray, tdipdip);
#endif
} TINKER_NAMESPACE_END

#endif
