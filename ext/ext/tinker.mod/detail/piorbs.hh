#ifndef TINKER_MOD_PIORBS_HH_
#define TINKER_MOD_PIORBS_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace piorbs {
extern int& norbit;
extern int& nconj;
extern int& reorbit;
extern int& nbpi;
extern int& ntpi;
extern int*& iorbit;
extern int*& iconj;
extern int*& kconj;
extern int*& piperp;
extern int*& ibpi;
extern int*& itpi;
extern double*& pbpl;
extern double*& pnpl;
extern int*& listpi;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(piorbs, norbit);
extern "C" int m_tinker_mod(piorbs, nconj);
extern "C" int m_tinker_mod(piorbs, reorbit);
extern "C" int m_tinker_mod(piorbs, nbpi);
extern "C" int m_tinker_mod(piorbs, ntpi);
extern "C" int* m_tinker_mod(piorbs, iorbit);
extern "C" int* m_tinker_mod(piorbs, iconj);
extern "C" int* m_tinker_mod(piorbs, kconj);
extern "C" int* m_tinker_mod(piorbs, piperp);
extern "C" int* m_tinker_mod(piorbs, ibpi);
extern "C" int* m_tinker_mod(piorbs, itpi);
extern "C" double* m_tinker_mod(piorbs, pbpl);
extern "C" double* m_tinker_mod(piorbs, pnpl);
extern "C" int* m_tinker_mod(piorbs, listpi);

int& norbit = m_tinker_mod(piorbs, norbit);
int& nconj = m_tinker_mod(piorbs, nconj);
int& reorbit = m_tinker_mod(piorbs, reorbit);
int& nbpi = m_tinker_mod(piorbs, nbpi);
int& ntpi = m_tinker_mod(piorbs, ntpi);
int*& iorbit = m_tinker_mod(piorbs, iorbit);
int*& iconj = m_tinker_mod(piorbs, iconj);
int*& kconj = m_tinker_mod(piorbs, kconj);
int*& piperp = m_tinker_mod(piorbs, piperp);
int*& ibpi = m_tinker_mod(piorbs, ibpi);
int*& itpi = m_tinker_mod(piorbs, itpi);
double*& pbpl = m_tinker_mod(piorbs, pbpl);
double*& pnpl = m_tinker_mod(piorbs, pnpl);
int*& listpi = m_tinker_mod(piorbs, listpi);
#endif

} TINKER_NAMESPACE_END

#endif
