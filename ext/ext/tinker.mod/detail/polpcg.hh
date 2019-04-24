#ifndef TINKER_MOD_POLPCG_HH_
#define TINKER_MOD_POLPCG_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace polpcg {
extern int*& mindex;
extern double& pcgpeek;
extern double*& minv;
extern int& pcgprec;
extern int& pcgguess;

#ifdef TINKER_MOD_CPP_
extern "C" int* m_tinker_mod(polpcg, mindex);
extern "C" double m_tinker_mod(polpcg, pcgpeek);
extern "C" double* m_tinker_mod(polpcg, minv);
extern "C" int m_tinker_mod(polpcg, pcgprec);
extern "C" int m_tinker_mod(polpcg, pcgguess);

int*& mindex = m_tinker_mod(polpcg, mindex);
double& pcgpeek = m_tinker_mod(polpcg, pcgpeek);
double*& minv = m_tinker_mod(polpcg, minv);
int& pcgprec = m_tinker_mod(polpcg, pcgprec);
int& pcgguess = m_tinker_mod(polpcg, pcgguess);
#endif

} TINKER_NAMESPACE_END

#endif
