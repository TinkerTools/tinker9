#ifndef TINKER_MOD_KREPL_HH_
#define TINKER_MOD_KREPL_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace krepl {
extern double*& prsiz;
extern double*& prdmp;
extern double*& prele;

#ifdef TINKER_MOD_CPP_
extern "C" double* m_tinker_mod(krepl, prsiz);
extern "C" double* m_tinker_mod(krepl, prdmp);
extern "C" double* m_tinker_mod(krepl, prele);

double*& prsiz = m_tinker_mod(krepl, prsiz);
double*& prdmp = m_tinker_mod(krepl, prdmp);
double*& prele = m_tinker_mod(krepl, prele);
#endif

} TINKER_NAMESPACE_END

#endif
