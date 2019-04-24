#ifndef TINKER_MOD_BITOR_HH_
#define TINKER_MOD_BITOR_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace bitor_ {
extern int& nbitor;
extern int*& ibitor;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(bitor, nbitor);
extern "C" int* m_tinker_mod(bitor, ibitor);

int& nbitor = m_tinker_mod(bitor, nbitor);
int*& ibitor = m_tinker_mod(bitor, ibitor);
#endif

} TINKER_NAMESPACE_END

#endif
