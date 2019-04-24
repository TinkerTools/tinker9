#ifndef TINKER_MOD_STRBND_HH_
#define TINKER_MOD_STRBND_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace strbnd {
extern int& nstrbnd;
extern int*& isb;
extern double*& sbk;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(strbnd, nstrbnd);
extern "C" int* m_tinker_mod(strbnd, isb);
extern "C" double* m_tinker_mod(strbnd, sbk);

int& nstrbnd = m_tinker_mod(strbnd, nstrbnd);
int*& isb = m_tinker_mod(strbnd, isb);
double*& sbk = m_tinker_mod(strbnd, sbk);
#endif

} TINKER_NAMESPACE_END

#endif
