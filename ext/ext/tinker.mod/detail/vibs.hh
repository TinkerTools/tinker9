#ifndef TINKER_MOD_VIBS_HH_
#define TINKER_MOD_VIBS_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace vibs {
extern double*& phi;
extern double*& phik;
extern double*& pwork;

#ifdef TINKER_MOD_CPP_
extern "C" double* m_tinker_mod(vibs, phi);
extern "C" double* m_tinker_mod(vibs, phik);
extern "C" double* m_tinker_mod(vibs, pwork);

double*& phi = m_tinker_mod(vibs, phi);
double*& phik = m_tinker_mod(vibs, phik);
double*& pwork = m_tinker_mod(vibs, pwork);
#endif

} TINKER_NAMESPACE_END

#endif
