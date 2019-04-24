#ifndef TINKER_MOD_KVDWS_HH_
#define TINKER_MOD_KVDWS_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kvdws {
extern double*& rad;
extern double*& eps;
extern double*& rad4;
extern double*& eps4;
extern double*& reduct;

#ifdef TINKER_MOD_CPP_
extern "C" double* m_tinker_mod(kvdws, rad);
extern "C" double* m_tinker_mod(kvdws, eps);
extern "C" double* m_tinker_mod(kvdws, rad4);
extern "C" double* m_tinker_mod(kvdws, eps4);
extern "C" double* m_tinker_mod(kvdws, reduct);

double*& rad = m_tinker_mod(kvdws, rad);
double*& eps = m_tinker_mod(kvdws, eps);
double*& rad4 = m_tinker_mod(kvdws, rad4);
double*& eps4 = m_tinker_mod(kvdws, eps4);
double*& reduct = m_tinker_mod(kvdws, reduct);
#endif

} TINKER_NAMESPACE_END

#endif
