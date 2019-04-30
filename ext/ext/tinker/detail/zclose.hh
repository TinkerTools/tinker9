#ifndef TINKER_MOD_ZCLOSE_HH_
#define TINKER_MOD_ZCLOSE_HH_

#include "util/macro.h"
#include "sizes.hh"

TINKER_NAMESPACE_BEGIN namespace zclose {
using namespace sizes;

extern int& nadd;
extern int& ndel;
extern int (&iadd)[maxatm][2];
extern int (&idel)[maxatm][2];

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(zclose, nadd);
extern "C" int TINKER_MOD(zclose, ndel);
extern "C" int TINKER_MOD(zclose, iadd)[maxatm][2];
extern "C" int TINKER_MOD(zclose, idel)[maxatm][2];

int& nadd = TINKER_MOD(zclose, nadd);
int& ndel = TINKER_MOD(zclose, ndel);
int (&iadd)[maxatm][2] = TINKER_MOD(zclose, iadd);
int (&idel)[maxatm][2] = TINKER_MOD(zclose, idel);
#endif
} TINKER_NAMESPACE_END

#endif
