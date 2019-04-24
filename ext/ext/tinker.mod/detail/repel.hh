#ifndef TINKER_MOD_REPEL_HH_
#define TINKER_MOD_REPEL_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace repel {
extern int& nrep;
extern double*& sizpr;
extern double*& dmppr;
extern double*& elepr;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(repel, nrep);
extern "C" double* m_tinker_mod(repel, sizpr);
extern "C" double* m_tinker_mod(repel, dmppr);
extern "C" double* m_tinker_mod(repel, elepr);

int& nrep = m_tinker_mod(repel, nrep);
double*& sizpr = m_tinker_mod(repel, sizpr);
double*& dmppr = m_tinker_mod(repel, dmppr);
double*& elepr = m_tinker_mod(repel, elepr);
#endif

} TINKER_NAMESPACE_END

#endif
