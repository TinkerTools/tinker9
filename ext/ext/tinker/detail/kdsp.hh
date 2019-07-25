#ifndef TINKER_MOD_KDSP_HH_
#define TINKER_MOD_KDSP_HH_

#include "util_macro.h"

TINKER_NAMESPACE_BEGIN namespace kdsp {
extern double*& dspsix;
extern double*& dspdmp;

#ifdef TINKER_MOD_CPP_
extern "C" double* TINKER_MOD(kdsp, dspsix);
extern "C" double* TINKER_MOD(kdsp, dspdmp);

double*& dspsix = TINKER_MOD(kdsp, dspsix);
double*& dspdmp = TINKER_MOD(kdsp, dspdmp);
#endif
} TINKER_NAMESPACE_END

#endif
