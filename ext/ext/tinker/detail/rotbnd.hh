#ifndef TINKER_MOD_ROTBND_HH_
#define TINKER_MOD_ROTBND_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace rotbnd {
extern int& nrot;
extern int*& rot;
extern int& use_short;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(rotbnd, nrot);
extern "C" int* TINKER_MOD(rotbnd, rot);
extern "C" int TINKER_MOD(rotbnd, use_short);

int& nrot = TINKER_MOD(rotbnd, nrot);
int*& rot = TINKER_MOD(rotbnd, rot);
int& use_short = TINKER_MOD(rotbnd, use_short);
#endif
} TINKER_NAMESPACE_END

#endif
