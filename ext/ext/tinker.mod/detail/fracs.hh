#ifndef TINKER_MOD_FRACS_HH_
#define TINKER_MOD_FRACS_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace fracs {
extern double*& xfrac;
extern double*& yfrac;
extern double*& zfrac;

#ifdef TINKER_MOD_CPP_
extern "C" double* m_tinker_mod(fracs, xfrac);
extern "C" double* m_tinker_mod(fracs, yfrac);
extern "C" double* m_tinker_mod(fracs, zfrac);

double*& xfrac = m_tinker_mod(fracs, xfrac);
double*& yfrac = m_tinker_mod(fracs, yfrac);
double*& zfrac = m_tinker_mod(fracs, zfrac);
#endif

} TINKER_NAMESPACE_END

#endif
