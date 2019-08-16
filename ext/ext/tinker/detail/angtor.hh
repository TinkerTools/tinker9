#ifndef TINKER_MOD_ANGTOR_HH_
#define TINKER_MOD_ANGTOR_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace angtor {
extern int& nangtor;
extern int*& iat;
extern double*& kant;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(angtor, nangtor);
extern "C" int* TINKER_MOD(angtor, iat);
extern "C" double* TINKER_MOD(angtor, kant);

int& nangtor = TINKER_MOD(angtor, nangtor);
int*& iat = TINKER_MOD(angtor, iat);
double*& kant = TINKER_MOD(angtor, kant);
#endif
} TINKER_NAMESPACE_END

#endif
