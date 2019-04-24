#ifndef TINKER_MOD_HESSN_HH_
#define TINKER_MOD_HESSN_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace hessn {
extern double*& hessx;
extern double*& hessy;
extern double*& hessz;

#ifdef TINKER_MOD_CPP_
extern "C" double* m_tinker_mod(hessn, hessx);
extern "C" double* m_tinker_mod(hessn, hessy);
extern "C" double* m_tinker_mod(hessn, hessz);

double*& hessx = m_tinker_mod(hessn, hessx);
double*& hessy = m_tinker_mod(hessn, hessy);
double*& hessz = m_tinker_mod(hessn, hessz);
#endif

} TINKER_NAMESPACE_END

#endif
