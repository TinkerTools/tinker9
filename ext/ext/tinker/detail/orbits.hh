#ifndef TINKER_MOD_ORBITS_HH_
#define TINKER_MOD_ORBITS_HH_

#include "util_macro.h"

TINKER_NAMESPACE_BEGIN namespace orbits {
extern double*& qorb;
extern double*& worb;
extern double*& emorb;

#ifdef TINKER_MOD_CPP_
extern "C" double* TINKER_MOD(orbits, qorb);
extern "C" double* TINKER_MOD(orbits, worb);
extern "C" double* TINKER_MOD(orbits, emorb);

double*& qorb = TINKER_MOD(orbits, qorb);
double*& worb = TINKER_MOD(orbits, worb);
double*& emorb = TINKER_MOD(orbits, emorb);
#endif
} TINKER_NAMESPACE_END

#endif
