#ifndef TINKER_MOD_RIGID_HH_
#define TINKER_MOD_RIGID_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace rigid {
extern double*& xrb;
extern double*& yrb;
extern double*& zrb;
extern double*& rbc;
extern int& use_rigid;

#ifdef TINKER_MOD_CPP_
extern "C" double* m_tinker_mod(rigid, xrb);
extern "C" double* m_tinker_mod(rigid, yrb);
extern "C" double* m_tinker_mod(rigid, zrb);
extern "C" double* m_tinker_mod(rigid, rbc);
extern "C" int m_tinker_mod(rigid, use_rigid);

double*& xrb = m_tinker_mod(rigid, xrb);
double*& yrb = m_tinker_mod(rigid, yrb);
double*& zrb = m_tinker_mod(rigid, zrb);
double*& rbc = m_tinker_mod(rigid, rbc);
int& use_rigid = m_tinker_mod(rigid, use_rigid);
#endif

} TINKER_NAMESPACE_END

#endif
