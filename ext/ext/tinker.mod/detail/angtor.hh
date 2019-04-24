#ifndef TINKER_MOD_ANGTOR_HH_
#define TINKER_MOD_ANGTOR_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace angtor {
extern int& nangtor;
extern int*& iat;
extern double*& kant;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(angtor, nangtor);
extern "C" int* m_tinker_mod(angtor, iat);
extern "C" double* m_tinker_mod(angtor, kant);

int& nangtor = m_tinker_mod(angtor, nangtor);
int*& iat = m_tinker_mod(angtor, iat);
double*& kant = m_tinker_mod(angtor, kant);
#endif

} TINKER_NAMESPACE_END

#endif
