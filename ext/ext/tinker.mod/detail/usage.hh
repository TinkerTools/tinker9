#ifndef TINKER_MOD_USAGE_HH_
#define TINKER_MOD_USAGE_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace usage {
extern int& nuse;
extern int*& iuse;
extern int*& use;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(usage, nuse);
extern "C" int* m_tinker_mod(usage, iuse);
extern "C" int* m_tinker_mod(usage, use);

int& nuse = m_tinker_mod(usage, nuse);
int*& iuse = m_tinker_mod(usage, iuse);
int*& use = m_tinker_mod(usage, use);
#endif

} TINKER_NAMESPACE_END

#endif
