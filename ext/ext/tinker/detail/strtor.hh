#ifndef TINKER_MOD_STRTOR_HH_
#define TINKER_MOD_STRTOR_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace strtor {
extern int& nstrtor;
extern int*& ist;
extern double*& kst;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(strtor, nstrtor);
extern "C" int* TINKER_MOD(strtor, ist);
extern "C" double* TINKER_MOD(strtor, kst);

int& nstrtor = TINKER_MOD(strtor, nstrtor);
int*& ist = TINKER_MOD(strtor, ist);
double*& kst = TINKER_MOD(strtor, kst);
#endif
} TINKER_NAMESPACE_END

#endif
