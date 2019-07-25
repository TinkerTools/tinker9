#ifndef TINKER_MOD_ATMLST_HH_
#define TINKER_MOD_ATMLST_HH_

#include "util_macro.h"

TINKER_NAMESPACE_BEGIN namespace atmlst {
extern int*& bndlist;
extern int*& anglist;

#ifdef TINKER_MOD_CPP_
extern "C" int* TINKER_MOD(atmlst, bndlist);
extern "C" int* TINKER_MOD(atmlst, anglist);

int*& bndlist = TINKER_MOD(atmlst, bndlist);
int*& anglist = TINKER_MOD(atmlst, anglist);
#endif
} TINKER_NAMESPACE_END

#endif
