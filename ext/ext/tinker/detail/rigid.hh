#ifndef TINKER_MOD_RIGID_HH_
#define TINKER_MOD_RIGID_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace rigid {
extern double*& xrb;
extern double*& yrb;
extern double*& zrb;
extern double*& rbc;
extern int& use_rigid;

#ifdef TINKER_MOD_CPP_
extern "C" double* TINKER_MOD(rigid, xrb);
extern "C" double* TINKER_MOD(rigid, yrb);
extern "C" double* TINKER_MOD(rigid, zrb);
extern "C" double* TINKER_MOD(rigid, rbc);
extern "C" int TINKER_MOD(rigid, use_rigid);

double*& xrb = TINKER_MOD(rigid, xrb);
double*& yrb = TINKER_MOD(rigid, yrb);
double*& zrb = TINKER_MOD(rigid, zrb);
double*& rbc = TINKER_MOD(rigid, rbc);
int& use_rigid = TINKER_MOD(rigid, use_rigid);
#endif
} TINKER_NAMESPACE_END

#endif
