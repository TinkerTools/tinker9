#ifndef TINKER_MOD_HESSN_HH_
#define TINKER_MOD_HESSN_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace hessn {
extern double*& hessx;
extern double*& hessy;
extern double*& hessz;

#ifdef TINKER_MOD_CPP_
extern "C" double* TINKER_MOD(hessn, hessx);
extern "C" double* TINKER_MOD(hessn, hessy);
extern "C" double* TINKER_MOD(hessn, hessz);

double*& hessx = TINKER_MOD(hessn, hessx);
double*& hessy = TINKER_MOD(hessn, hessy);
double*& hessz = TINKER_MOD(hessn, hessz);
#endif
} TINKER_NAMESPACE_END

#endif
