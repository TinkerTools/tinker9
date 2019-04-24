#ifndef TINKER_MOD_PISTUF_HH_
#define TINKER_MOD_PISTUF_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace pistuf {
extern double*& bkpi;
extern double*& blpi;
extern double*& kslope;
extern double*& lslope;
extern double*& torsp2;

#ifdef TINKER_MOD_CPP_
extern "C" double* m_tinker_mod(pistuf, bkpi);
extern "C" double* m_tinker_mod(pistuf, blpi);
extern "C" double* m_tinker_mod(pistuf, kslope);
extern "C" double* m_tinker_mod(pistuf, lslope);
extern "C" double* m_tinker_mod(pistuf, torsp2);

double*& bkpi = m_tinker_mod(pistuf, bkpi);
double*& blpi = m_tinker_mod(pistuf, blpi);
double*& kslope = m_tinker_mod(pistuf, kslope);
double*& lslope = m_tinker_mod(pistuf, lslope);
double*& torsp2 = m_tinker_mod(pistuf, torsp2);
#endif

} TINKER_NAMESPACE_END

#endif
