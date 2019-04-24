#ifndef TINKER_MOD_TORTOR_HH_
#define TINKER_MOD_TORTOR_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace tortor {
extern int& ntortor;
extern int*& itt;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(tortor, ntortor);
extern "C" int* m_tinker_mod(tortor, itt);

int& ntortor = m_tinker_mod(tortor, ntortor);
int*& itt = m_tinker_mod(tortor, itt);
#endif

} TINKER_NAMESPACE_END

#endif
