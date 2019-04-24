#ifndef TINKER_MOD_POLTCG_HH_
#define TINKER_MOD_POLTCG_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace poltcg {
extern int& tcgorder;
extern int& tcgnab;
extern double& tcgpeek;
extern double*& uad;
extern double*& uap;
extern double*& ubd;
extern double*& ubp;
extern int& tcgguess;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(poltcg, tcgorder);
extern "C" int m_tinker_mod(poltcg, tcgnab);
extern "C" double m_tinker_mod(poltcg, tcgpeek);
extern "C" double* m_tinker_mod(poltcg, uad);
extern "C" double* m_tinker_mod(poltcg, uap);
extern "C" double* m_tinker_mod(poltcg, ubd);
extern "C" double* m_tinker_mod(poltcg, ubp);
extern "C" int m_tinker_mod(poltcg, tcgguess);

int& tcgorder = m_tinker_mod(poltcg, tcgorder);
int& tcgnab = m_tinker_mod(poltcg, tcgnab);
double& tcgpeek = m_tinker_mod(poltcg, tcgpeek);
double*& uad = m_tinker_mod(poltcg, uad);
double*& uap = m_tinker_mod(poltcg, uap);
double*& ubd = m_tinker_mod(poltcg, ubd);
double*& ubp = m_tinker_mod(poltcg, ubp);
int& tcgguess = m_tinker_mod(poltcg, tcgguess);
#endif

} TINKER_NAMESPACE_END

#endif
