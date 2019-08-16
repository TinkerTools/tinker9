#ifndef TINKER_MOD_IOUNIT_HH_
#define TINKER_MOD_IOUNIT_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace iounit {
extern int& input;
extern int& iout;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(iounit, input);
extern "C" int TINKER_MOD(iounit, iout);

int& input = TINKER_MOD(iounit, input);
int& iout = TINKER_MOD(iounit, iout);
#endif
} TINKER_NAMESPACE_END

#endif
