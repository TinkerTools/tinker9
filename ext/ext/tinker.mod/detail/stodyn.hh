#ifndef TINKER_MOD_STODYN_HH_
#define TINKER_MOD_STODYN_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace stodyn {
extern double& friction;
extern double*& fgamma;
extern int& use_sdarea;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(stodyn, friction);
extern "C" double* m_tinker_mod(stodyn, fgamma);
extern "C" int m_tinker_mod(stodyn, use_sdarea);

double& friction = m_tinker_mod(stodyn, friction);
double*& fgamma = m_tinker_mod(stodyn, fgamma);
int& use_sdarea = m_tinker_mod(stodyn, use_sdarea);
#endif

} TINKER_NAMESPACE_END

#endif
