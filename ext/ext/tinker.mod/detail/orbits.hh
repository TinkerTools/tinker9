#ifndef TINKER_MOD_ORBITS_HH_
#define TINKER_MOD_ORBITS_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace orbits {
extern double*& qorb;
extern double*& worb;
extern double*& emorb;

#ifdef TINKER_MOD_CPP_
extern "C" double* m_tinker_mod(orbits, qorb);
extern "C" double* m_tinker_mod(orbits, worb);
extern "C" double* m_tinker_mod(orbits, emorb);

double*& qorb = m_tinker_mod(orbits, qorb);
double*& worb = m_tinker_mod(orbits, worb);
double*& emorb = m_tinker_mod(orbits, emorb);
#endif

} TINKER_NAMESPACE_END

#endif
