#ifndef TINKER_MOD_STRBND_HH_
#define TINKER_MOD_STRBND_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace strbnd {
extern int& nstrbnd;
extern int*& isb;
extern double*& sbk;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(strbnd, nstrbnd);
extern "C" int* TINKER_MOD(strbnd, isb);
extern "C" double* TINKER_MOD(strbnd, sbk);

int& nstrbnd = TINKER_MOD(strbnd, nstrbnd);
int*& isb = TINKER_MOD(strbnd, isb);
double*& sbk = TINKER_MOD(strbnd, sbk);
#endif
} TINKER_NAMESPACE_END

#endif
