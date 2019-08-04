#ifndef TINKER_MOD_BNDSTR_HH_
#define TINKER_MOD_BNDSTR_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace bndstr {
extern int& nbond;
extern int*& ibnd;
extern double*& bk;
extern double*& bl;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(bndstr, nbond);
extern "C" int* TINKER_MOD(bndstr, ibnd);
extern "C" double* TINKER_MOD(bndstr, bk);
extern "C" double* TINKER_MOD(bndstr, bl);

int& nbond = TINKER_MOD(bndstr, nbond);
int*& ibnd = TINKER_MOD(bndstr, ibnd);
double*& bk = TINKER_MOD(bndstr, bk);
double*& bl = TINKER_MOD(bndstr, bl);
#endif
} TINKER_NAMESPACE_END

#endif
