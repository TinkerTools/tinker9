#ifndef TINKER_MOD_MOLDYN_HH_
#define TINKER_MOD_MOLDYN_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace moldyn {
extern double*& v;
extern double*& a;
extern double*& aalt;

#ifdef TINKER_MOD_CPP_
extern "C" double* m_tinker_mod(moldyn, v);
extern "C" double* m_tinker_mod(moldyn, a);
extern "C" double* m_tinker_mod(moldyn, aalt);

double*& v = m_tinker_mod(moldyn, v);
double*& a = m_tinker_mod(moldyn, a);
double*& aalt = m_tinker_mod(moldyn, aalt);
#endif

} TINKER_NAMESPACE_END

#endif
