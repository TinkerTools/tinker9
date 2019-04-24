#ifndef TINKER_MOD_MINIMA_HH_
#define TINKER_MOD_MINIMA_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace minima {
extern int& maxiter;
extern int& nextiter;
extern double& fctmin;
extern double& hguess;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(minima, maxiter);
extern "C" int m_tinker_mod(minima, nextiter);
extern "C" double m_tinker_mod(minima, fctmin);
extern "C" double m_tinker_mod(minima, hguess);

int& maxiter = m_tinker_mod(minima, maxiter);
int& nextiter = m_tinker_mod(minima, nextiter);
double& fctmin = m_tinker_mod(minima, fctmin);
double& hguess = m_tinker_mod(minima, hguess);
#endif

} TINKER_NAMESPACE_END

#endif
