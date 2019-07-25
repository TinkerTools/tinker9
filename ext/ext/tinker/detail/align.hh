#ifndef TINKER_MOD_ALIGN_HH_
#define TINKER_MOD_ALIGN_HH_

#include "util_macro.h"

TINKER_NAMESPACE_BEGIN namespace align {
extern int& nfit;
extern int*& ifit;
extern double*& wfit;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(align, nfit);
extern "C" int* TINKER_MOD(align, ifit);
extern "C" double* TINKER_MOD(align, wfit);

int& nfit = TINKER_MOD(align, nfit);
int*& ifit = TINKER_MOD(align, ifit);
double*& wfit = TINKER_MOD(align, wfit);
#endif
} TINKER_NAMESPACE_END

#endif
