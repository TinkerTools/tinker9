#ifndef TINKER_MOD_OPBEND_HH_
#define TINKER_MOD_OPBEND_HH_

#include "util_macro.h"

TINKER_NAMESPACE_BEGIN namespace opbend {
extern int& nopbend;
extern int*& iopb;
extern double*& opbk;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(opbend, nopbend);
extern "C" int* TINKER_MOD(opbend, iopb);
extern "C" double* TINKER_MOD(opbend, opbk);

int& nopbend = TINKER_MOD(opbend, nopbend);
int*& iopb = TINKER_MOD(opbend, iopb);
double*& opbk = TINKER_MOD(opbend, opbk);
#endif
} TINKER_NAMESPACE_END

#endif
