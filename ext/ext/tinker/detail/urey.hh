#ifndef TINKER_MOD_UREY_HH_
#define TINKER_MOD_UREY_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace urey {
extern int& nurey;
extern int*& iury;
extern double*& uk;
extern double*& ul;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(urey, nurey);
extern "C" int* TINKER_MOD(urey, iury);
extern "C" double* TINKER_MOD(urey, uk);
extern "C" double* TINKER_MOD(urey, ul);

int& nurey = TINKER_MOD(urey, nurey);
int*& iury = TINKER_MOD(urey, iury);
double*& uk = TINKER_MOD(urey, uk);
double*& ul = TINKER_MOD(urey, ul);
#endif

} TINKER_NAMESPACE_END

#endif
