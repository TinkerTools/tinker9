#ifndef TINKER_MOD_VIBS_HH_
#define TINKER_MOD_VIBS_HH_

#include "util_macro.h"

TINKER_NAMESPACE_BEGIN namespace vibs {
extern double*& phi;
extern double*& phik;
extern double*& pwork;

#ifdef TINKER_MOD_CPP_
extern "C" double* TINKER_MOD(vibs, phi);
extern "C" double* TINKER_MOD(vibs, phik);
extern "C" double* TINKER_MOD(vibs, pwork);

double*& phi = TINKER_MOD(vibs, phi);
double*& phik = TINKER_MOD(vibs, phik);
double*& pwork = TINKER_MOD(vibs, pwork);
#endif
} TINKER_NAMESPACE_END

#endif
