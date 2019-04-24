#ifndef TINKER_MOD_KDSP_HH_
#define TINKER_MOD_KDSP_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kdsp {
extern double*& dspsix;
extern double*& dspdmp;

#ifdef TINKER_MOD_CPP_
extern "C" double* m_tinker_mod(kdsp, dspsix);
extern "C" double* m_tinker_mod(kdsp, dspdmp);

double*& dspsix = m_tinker_mod(kdsp, dspsix);
double*& dspdmp = m_tinker_mod(kdsp, dspdmp);
#endif

} TINKER_NAMESPACE_END

#endif
