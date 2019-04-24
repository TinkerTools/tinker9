#ifndef TINKER_MOD_IOUNIT_HH_
#define TINKER_MOD_IOUNIT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace iounit {
extern int& input;
extern int& iout;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(iounit, input);
extern "C" int m_tinker_mod(iounit, iout);

int& input = m_tinker_mod(iounit, input);
int& iout = m_tinker_mod(iounit, iout);
#endif

} TINKER_NAMESPACE_END

#endif
