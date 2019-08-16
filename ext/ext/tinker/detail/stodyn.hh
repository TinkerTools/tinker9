#ifndef TINKER_MOD_STODYN_HH_
#define TINKER_MOD_STODYN_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace stodyn {
extern double& friction;
extern double*& fgamma;
extern int& use_sdarea;

#ifdef TINKER_MOD_CPP_
extern "C" double TINKER_MOD(stodyn, friction);
extern "C" double* TINKER_MOD(stodyn, fgamma);
extern "C" int TINKER_MOD(stodyn, use_sdarea);

double& friction = TINKER_MOD(stodyn, friction);
double*& fgamma = TINKER_MOD(stodyn, fgamma);
int& use_sdarea = TINKER_MOD(stodyn, use_sdarea);
#endif
} TINKER_NAMESPACE_END

#endif
