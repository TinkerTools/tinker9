#ifndef TINKER_MOD_UREY_HH_
#define TINKER_MOD_UREY_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace urey {
extern int& nurey;
extern int*& iury;
extern double*& uk;
extern double*& ul;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(urey, nurey);
extern "C" int* m_tinker_mod(urey, iury);
extern "C" double* m_tinker_mod(urey, uk);
extern "C" double* m_tinker_mod(urey, ul);

int& nurey = m_tinker_mod(urey, nurey);
int*& iury = m_tinker_mod(urey, iury);
double*& uk = m_tinker_mod(urey, uk);
double*& ul = m_tinker_mod(urey, ul);
#endif

} TINKER_NAMESPACE_END

#endif
