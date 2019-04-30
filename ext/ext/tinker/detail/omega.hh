#ifndef TINKER_MOD_OMEGA_HH_
#define TINKER_MOD_OMEGA_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace omega {
extern int& nomega;
extern int*& iomega;
extern int*& zline;
extern double*& dihed;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(omega, nomega);
extern "C" int* TINKER_MOD(omega, iomega);
extern "C" int* TINKER_MOD(omega, zline);
extern "C" double* TINKER_MOD(omega, dihed);

int& nomega = TINKER_MOD(omega, nomega);
int*& iomega = TINKER_MOD(omega, iomega);
int*& zline = TINKER_MOD(omega, zline);
double*& dihed = TINKER_MOD(omega, dihed);
#endif
} TINKER_NAMESPACE_END

#endif
