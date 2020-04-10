#ifndef TINKER_MOD_VIBS_HH_
#define TINKER_MOD_VIBS_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace vibs {
extern double*& rho;
extern double*& rhok;
extern double*& rwork;

#ifdef TINKER_MOD_CPP_
extern "C" double* TINKER_MOD(vibs, rho);
extern "C" double* TINKER_MOD(vibs, rhok);
extern "C" double* TINKER_MOD(vibs, rwork);

double*& rho = TINKER_MOD(vibs, rho);
double*& rhok = TINKER_MOD(vibs, rhok);
double*& rwork = TINKER_MOD(vibs, rwork);
#endif
} TINKER_NAMESPACE_END

#endif
