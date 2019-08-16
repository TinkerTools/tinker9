#ifndef TINKER_MOD_FRACS_HH_
#define TINKER_MOD_FRACS_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace fracs {
extern double*& xfrac;
extern double*& yfrac;
extern double*& zfrac;

#ifdef TINKER_MOD_CPP_
extern "C" double* TINKER_MOD(fracs, xfrac);
extern "C" double* TINKER_MOD(fracs, yfrac);
extern "C" double* TINKER_MOD(fracs, zfrac);

double*& xfrac = TINKER_MOD(fracs, xfrac);
double*& yfrac = TINKER_MOD(fracs, yfrac);
double*& zfrac = TINKER_MOD(fracs, zfrac);
#endif
} TINKER_NAMESPACE_END

#endif
