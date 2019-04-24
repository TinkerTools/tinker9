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
extern "C" int m_tinker_mod(zclose, nadd);
extern "C" int m_tinker_mod(zclose, ndel);
extern "C" int m_tinker_mod(zclose, iadd)[maxatm][2];
extern "C" int m_tinker_mod(zclose, idel)[maxatm][2];

int& nadd = m_tinker_mod(zclose, nadd);
int& ndel = m_tinker_mod(zclose, ndel);
int (&iadd)[maxatm][2] = m_tinker_mod(zclose, iadd);
int (&idel)[maxatm][2] = m_tinker_mod(zclose, idel);
#endif

} TINKER_NAMESPACE_END

#endif
