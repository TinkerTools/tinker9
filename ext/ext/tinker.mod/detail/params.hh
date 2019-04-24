#ifndef TINKER_MOD_PARAMS_HH_
#define TINKER_MOD_PARAMS_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace params {
const int maxprm = 25000;
extern int& nprm;
extern char (&prmline)[maxprm][240];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(params, nprm);
extern "C" char m_tinker_mod(params, prmline)[maxprm][240];

int& nprm = m_tinker_mod(params, nprm);
char (&prmline)[maxprm][240] = m_tinker_mod(params, prmline);
#endif

} TINKER_NAMESPACE_END

#endif
