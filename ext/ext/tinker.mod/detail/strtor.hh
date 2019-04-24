#ifndef TINKER_MOD_STRTOR_HH_
#define TINKER_MOD_STRTOR_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace strtor {
extern int& nstrtor;
extern int*& ist;
extern double*& kst;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(strtor, nstrtor);
extern "C" int* m_tinker_mod(strtor, ist);
extern "C" double* m_tinker_mod(strtor, kst);

int& nstrtor = m_tinker_mod(strtor, nstrtor);
int*& ist = m_tinker_mod(strtor, ist);
double*& kst = m_tinker_mod(strtor, kst);
#endif

} TINKER_NAMESPACE_END

#endif
