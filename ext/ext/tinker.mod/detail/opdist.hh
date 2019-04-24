#ifndef TINKER_MOD_OPDIST_HH_
#define TINKER_MOD_OPDIST_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace opdist {
extern int& nopdist;
extern int*& iopd;
extern double*& opdk;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(opdist, nopdist);
extern "C" int* TINKER_MOD(opdist, iopd);
extern "C" double* TINKER_MOD(opdist, opdk);

int& nopdist = TINKER_MOD(opdist, nopdist);
int*& iopd = TINKER_MOD(opdist, iopd);
double*& opdk = TINKER_MOD(opdist, opdk);
#endif

} TINKER_NAMESPACE_END

#endif
