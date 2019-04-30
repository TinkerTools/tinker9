#ifndef TINKER_MOD_MINIMA_HH_
#define TINKER_MOD_MINIMA_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace minima {
extern int& maxiter;
extern int& nextiter;
extern double& fctmin;
extern double& hguess;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(minima, maxiter);
extern "C" int TINKER_MOD(minima, nextiter);
extern "C" double TINKER_MOD(minima, fctmin);
extern "C" double TINKER_MOD(minima, hguess);

int& maxiter = TINKER_MOD(minima, maxiter);
int& nextiter = TINKER_MOD(minima, nextiter);
double& fctmin = TINKER_MOD(minima, fctmin);
double& hguess = TINKER_MOD(minima, hguess);
#endif
} TINKER_NAMESPACE_END

#endif
