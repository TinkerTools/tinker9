#pragma once

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace bitor_ {
extern int& nbitor;
extern int*& ibitor;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(bitor, nbitor);
extern "C" int* TINKER_MOD(bitor, ibitor);

int& nbitor = TINKER_MOD(bitor, nbitor);
int*& ibitor = TINKER_MOD(bitor, ibitor);
#endif
} TINKER_NAMESPACE_END
