#ifndef TINKER_MOD_USAGE_HH_
#define TINKER_MOD_USAGE_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace usage {
extern int& nuse;
extern int*& iuse;
extern int*& use;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(usage, nuse);
extern "C" int* TINKER_MOD(usage, iuse);
extern "C" int* TINKER_MOD(usage, use);

int& nuse = TINKER_MOD(usage, nuse);
int*& iuse = TINKER_MOD(usage, iuse);
int*& use = TINKER_MOD(usage, use);
#endif

} TINKER_NAMESPACE_END

#endif
