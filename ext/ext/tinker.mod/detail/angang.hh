#ifndef TINKER_MOD_ANGANG_HH_
#define TINKER_MOD_ANGANG_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace angang {
extern int& nangang;
extern int*& iaa;
extern double*& kaa;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(angang, nangang);
extern "C" int* m_tinker_mod(angang, iaa);
extern "C" double* m_tinker_mod(angang, kaa);

int& nangang = m_tinker_mod(angang, nangang);
int*& iaa = m_tinker_mod(angang, iaa);
double*& kaa = m_tinker_mod(angang, kaa);
#endif

} TINKER_NAMESPACE_END

#endif
