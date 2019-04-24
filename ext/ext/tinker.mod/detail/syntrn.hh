#ifndef TINKER_MOD_SYNTRN_HH_
#define TINKER_MOD_SYNTRN_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace syntrn {
extern double& tpath;
extern double& ppath;
extern double*& xmin1;
extern double*& xmin2;
extern double*& xm;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(syntrn, tpath);
extern "C" double m_tinker_mod(syntrn, ppath);
extern "C" double* m_tinker_mod(syntrn, xmin1);
extern "C" double* m_tinker_mod(syntrn, xmin2);
extern "C" double* m_tinker_mod(syntrn, xm);

double& tpath = m_tinker_mod(syntrn, tpath);
double& ppath = m_tinker_mod(syntrn, ppath);
double*& xmin1 = m_tinker_mod(syntrn, xmin1);
double*& xmin2 = m_tinker_mod(syntrn, xmin2);
double*& xm = m_tinker_mod(syntrn, xm);
#endif

} TINKER_NAMESPACE_END

#endif
