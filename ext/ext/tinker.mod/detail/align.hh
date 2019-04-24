#ifndef TINKER_MOD_ALIGN_HH_
#define TINKER_MOD_ALIGN_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace align {
extern int& nfit;
extern int*& ifit;
extern double*& wfit;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(align, nfit);
extern "C" int* m_tinker_mod(align, ifit);
extern "C" double* m_tinker_mod(align, wfit);

int& nfit = m_tinker_mod(align, nfit);
int*& ifit = m_tinker_mod(align, ifit);
double*& wfit = m_tinker_mod(align, wfit);
#endif

} TINKER_NAMESPACE_END

#endif
