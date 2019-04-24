#ifndef TINKER_MOD_DISP_HH_
#define TINKER_MOD_DISP_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace disp {
extern int& ndisp;
extern int*& idisp;
extern double& csixpr;
extern double*& csix;
extern double*& adisp;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(disp, ndisp);
extern "C" int* TINKER_MOD(disp, idisp);
extern "C" double TINKER_MOD(disp, csixpr);
extern "C" double* TINKER_MOD(disp, csix);
extern "C" double* TINKER_MOD(disp, adisp);

int& ndisp = TINKER_MOD(disp, ndisp);
int*& idisp = TINKER_MOD(disp, idisp);
double& csixpr = TINKER_MOD(disp, csixpr);
double*& csix = TINKER_MOD(disp, csix);
double*& adisp = TINKER_MOD(disp, adisp);
#endif

} TINKER_NAMESPACE_END

#endif
