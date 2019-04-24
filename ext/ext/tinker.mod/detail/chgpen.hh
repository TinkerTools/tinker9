#ifndef TINKER_MOD_CHGPEN_HH_
#define TINKER_MOD_CHGPEN_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace chgpen {
extern int& ncp;
extern double*& pcore;
extern double*& pval;
extern double*& palpha;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(chgpen, ncp);
extern "C" double* m_tinker_mod(chgpen, pcore);
extern "C" double* m_tinker_mod(chgpen, pval);
extern "C" double* m_tinker_mod(chgpen, palpha);

int& ncp = m_tinker_mod(chgpen, ncp);
double*& pcore = m_tinker_mod(chgpen, pcore);
double*& pval = m_tinker_mod(chgpen, pval);
double*& palpha = m_tinker_mod(chgpen, palpha);
#endif

} TINKER_NAMESPACE_END

#endif
