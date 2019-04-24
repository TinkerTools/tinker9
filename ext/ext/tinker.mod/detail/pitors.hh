#ifndef TINKER_MOD_PITORS_HH_
#define TINKER_MOD_PITORS_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace pitors {
extern int& npitors;
extern int*& ipit;
extern double*& kpit;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(pitors, npitors);
extern "C" int* m_tinker_mod(pitors, ipit);
extern "C" double* m_tinker_mod(pitors, kpit);

int& npitors = m_tinker_mod(pitors, npitors);
int*& ipit = m_tinker_mod(pitors, ipit);
double*& kpit = m_tinker_mod(pitors, kpit);
#endif

} TINKER_NAMESPACE_END

#endif
