#ifndef TINKER_MOD_OMEGA_HH_
#define TINKER_MOD_OMEGA_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace omega {
extern int& nomega;
extern int*& iomega;
extern int*& zline;
extern double*& dihed;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(omega, nomega);
extern "C" int* m_tinker_mod(omega, iomega);
extern "C" int* m_tinker_mod(omega, zline);
extern "C" double* m_tinker_mod(omega, dihed);

int& nomega = m_tinker_mod(omega, nomega);
int*& iomega = m_tinker_mod(omega, iomega);
int*& zline = m_tinker_mod(omega, zline);
double*& dihed = m_tinker_mod(omega, dihed);
#endif

} TINKER_NAMESPACE_END

#endif
