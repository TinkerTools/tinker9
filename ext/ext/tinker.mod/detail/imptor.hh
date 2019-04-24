#ifndef TINKER_MOD_IMPTOR_HH_
#define TINKER_MOD_IMPTOR_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace imptor {
extern int& nitors;
extern int*& iitors;
extern double*& itors1;
extern double*& itors2;
extern double*& itors3;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(imptor, nitors);
extern "C" int* m_tinker_mod(imptor, iitors);
extern "C" double* m_tinker_mod(imptor, itors1);
extern "C" double* m_tinker_mod(imptor, itors2);
extern "C" double* m_tinker_mod(imptor, itors3);

int& nitors = m_tinker_mod(imptor, nitors);
int*& iitors = m_tinker_mod(imptor, iitors);
double*& itors1 = m_tinker_mod(imptor, itors1);
double*& itors2 = m_tinker_mod(imptor, itors2);
double*& itors3 = m_tinker_mod(imptor, itors3);
#endif

} TINKER_NAMESPACE_END

#endif
