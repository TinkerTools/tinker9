#ifndef TINKER_MOD_MOLDYN_HH_
#define TINKER_MOD_MOLDYN_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace moldyn {
extern double*& v;
extern double*& a;
extern double*& aalt;

#ifdef TINKER_MOD_CPP_
extern "C" double* TINKER_MOD(moldyn, v);
extern "C" double* TINKER_MOD(moldyn, a);
extern "C" double* TINKER_MOD(moldyn, aalt);

double*& v = TINKER_MOD(moldyn, v);
double*& a = TINKER_MOD(moldyn, a);
double*& aalt = TINKER_MOD(moldyn, aalt);
#endif

} TINKER_NAMESPACE_END

#endif
