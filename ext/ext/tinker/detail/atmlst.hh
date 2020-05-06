#pragma once

#include "macro.h"

namespace tinker { namespace atmlst {
extern int*& bndlist;
extern int*& anglist;

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" int* TINKER_MOD(atmlst, bndlist);
extern "C" int* TINKER_MOD(atmlst, anglist);

int*& bndlist = TINKER_MOD(atmlst, bndlist);
int*& anglist = TINKER_MOD(atmlst, anglist);
#endif
} }
