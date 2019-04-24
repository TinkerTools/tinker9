#ifndef TINKER_MOD_KCPEN_HH_
#define TINKER_MOD_KCPEN_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kcpen {
extern double*& cpele;
extern double*& cpalp;

#ifdef TINKER_MOD_CPP_
extern "C" double* m_tinker_mod(kcpen, cpele);
extern "C" double* m_tinker_mod(kcpen, cpalp);

double*& cpele = m_tinker_mod(kcpen, cpele);
double*& cpalp = m_tinker_mod(kcpen, cpalp);
#endif

} TINKER_NAMESPACE_END

#endif
