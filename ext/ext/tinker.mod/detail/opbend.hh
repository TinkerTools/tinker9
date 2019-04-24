#ifndef TINKER_MOD_OPBEND_HH_
#define TINKER_MOD_OPBEND_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace opbend {
extern int& nopbend;
extern int*& iopb;
extern double*& opbk;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(opbend, nopbend);
extern "C" int* m_tinker_mod(opbend, iopb);
extern "C" double* m_tinker_mod(opbend, opbk);

int& nopbend = m_tinker_mod(opbend, nopbend);
int*& iopb = m_tinker_mod(opbend, iopb);
double*& opbk = m_tinker_mod(opbend, opbk);
#endif

} TINKER_NAMESPACE_END

#endif
