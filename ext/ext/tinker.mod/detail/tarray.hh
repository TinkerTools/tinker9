#ifndef TINKER_MOD_TARRAY_HH_
#define TINKER_MOD_TARRAY_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace tarray {
extern int& ntpair;
extern int*& tindex;
extern double*& tdipdip;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(tarray, ntpair);
extern "C" int* m_tinker_mod(tarray, tindex);
extern "C" double* m_tinker_mod(tarray, tdipdip);

int& ntpair = m_tinker_mod(tarray, ntpair);
int*& tindex = m_tinker_mod(tarray, tindex);
double*& tdipdip = m_tinker_mod(tarray, tdipdip);
#endif

} TINKER_NAMESPACE_END

#endif
