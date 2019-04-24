#ifndef TINKER_MOD_CHARGE_HH_
#define TINKER_MOD_CHARGE_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace charge {
extern int& nion;
extern int*& iion;
extern int*& jion;
extern int*& kion;
extern double*& pchg;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(charge, nion);
extern "C" int* m_tinker_mod(charge, iion);
extern "C" int* m_tinker_mod(charge, jion);
extern "C" int* m_tinker_mod(charge, kion);
extern "C" double* m_tinker_mod(charge, pchg);

int& nion = m_tinker_mod(charge, nion);
int*& iion = m_tinker_mod(charge, iion);
int*& jion = m_tinker_mod(charge, jion);
int*& kion = m_tinker_mod(charge, kion);
double*& pchg = m_tinker_mod(charge, pchg);
#endif

} TINKER_NAMESPACE_END

#endif
