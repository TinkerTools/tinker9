#ifndef TINKER_MOD_ANGBND_HH_
#define TINKER_MOD_ANGBND_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace angbnd {
extern int& nangle;
extern int*& iang;
extern double*& ak;
extern double*& anat;
extern double*& afld;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(angbnd, nangle);
extern "C" int* TINKER_MOD(angbnd, iang);
extern "C" double* TINKER_MOD(angbnd, ak);
extern "C" double* TINKER_MOD(angbnd, anat);
extern "C" double* TINKER_MOD(angbnd, afld);

int& nangle = TINKER_MOD(angbnd, nangle);
int*& iang = TINKER_MOD(angbnd, iang);
double*& ak = TINKER_MOD(angbnd, ak);
double*& anat = TINKER_MOD(angbnd, anat);
double*& afld = TINKER_MOD(angbnd, afld);
#endif

} TINKER_NAMESPACE_END

#endif
