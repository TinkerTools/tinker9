#ifndef TINKER_MOD_KREPL_HH_
#define TINKER_MOD_KREPL_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace krepl {
extern double*& prsiz;
extern double*& prdmp;
extern double*& prele;

#ifdef TINKER_MOD_CPP_
extern "C" double* TINKER_MOD(krepl, prsiz);
extern "C" double* TINKER_MOD(krepl, prdmp);
extern "C" double* TINKER_MOD(krepl, prele);

double*& prsiz = TINKER_MOD(krepl, prsiz);
double*& prdmp = TINKER_MOD(krepl, prdmp);
double*& prele = TINKER_MOD(krepl, prele);
#endif
} TINKER_NAMESPACE_END

#endif
