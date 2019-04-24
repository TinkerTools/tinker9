#ifndef TINKER_MOD_IMPROP_HH_
#define TINKER_MOD_IMPROP_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace improp {
extern int& niprop;
extern int*& iiprop;
extern double*& kprop;
extern double*& vprop;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(improp, niprop);
extern "C" int* m_tinker_mod(improp, iiprop);
extern "C" double* m_tinker_mod(improp, kprop);
extern "C" double* m_tinker_mod(improp, vprop);

int& niprop = m_tinker_mod(improp, niprop);
int*& iiprop = m_tinker_mod(improp, iiprop);
double*& kprop = m_tinker_mod(improp, kprop);
double*& vprop = m_tinker_mod(improp, vprop);
#endif

} TINKER_NAMESPACE_END

#endif
