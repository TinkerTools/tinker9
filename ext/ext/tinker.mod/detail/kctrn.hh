#ifndef TINKER_MOD_KCTRN_HH_
#define TINKER_MOD_KCTRN_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kctrn {
extern double*& ctchg;
extern double*& ctdmp;

#ifdef TINKER_MOD_CPP_
extern "C" double* m_tinker_mod(kctrn, ctchg);
extern "C" double* m_tinker_mod(kctrn, ctdmp);

double*& ctchg = m_tinker_mod(kctrn, ctchg);
double*& ctdmp = m_tinker_mod(kctrn, ctdmp);
#endif

} TINKER_NAMESPACE_END

#endif
