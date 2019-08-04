#ifndef TINKER_MOD_KCPEN_HH_
#define TINKER_MOD_KCPEN_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace kcpen {
extern double*& cpele;
extern double*& cpalp;

#ifdef TINKER_MOD_CPP_
extern "C" double* TINKER_MOD(kcpen, cpele);
extern "C" double* TINKER_MOD(kcpen, cpalp);

double*& cpele = TINKER_MOD(kcpen, cpele);
double*& cpalp = TINKER_MOD(kcpen, cpalp);
#endif
} TINKER_NAMESPACE_END

#endif
