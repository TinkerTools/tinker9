#ifndef TINKER_MOD_PARAMS_HH_
#define TINKER_MOD_PARAMS_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace params {
const int maxprm = 25000;
extern int& nprm;
extern char (&prmline)[maxprm][240];

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(params, nprm);
extern "C" char TINKER_MOD(params, prmline)[maxprm][240];

int& nprm = TINKER_MOD(params, nprm);
char (&prmline)[maxprm][240] = TINKER_MOD(params, prmline);
#endif
} TINKER_NAMESPACE_END

#endif
