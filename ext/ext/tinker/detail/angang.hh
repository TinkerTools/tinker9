#ifndef TINKER_MOD_ANGANG_HH_
#define TINKER_MOD_ANGANG_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace angang {
extern int& nangang;
extern int*& iaa;
extern double*& kaa;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(angang, nangang);
extern "C" int* TINKER_MOD(angang, iaa);
extern "C" double* TINKER_MOD(angang, kaa);

int& nangang = TINKER_MOD(angang, nangang);
int*& iaa = TINKER_MOD(angang, iaa);
double*& kaa = TINKER_MOD(angang, kaa);
#endif
} TINKER_NAMESPACE_END

#endif
