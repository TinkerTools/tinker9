#ifndef TINKER_MOD_ROTBND_HH_
#define TINKER_MOD_ROTBND_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace rotbnd {
extern int& nrot;
extern int*& rot;
extern int& use_short;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(rotbnd, nrot);
extern "C" int* m_tinker_mod(rotbnd, rot);
extern "C" int m_tinker_mod(rotbnd, use_short);

int& nrot = m_tinker_mod(rotbnd, nrot);
int*& rot = m_tinker_mod(rotbnd, rot);
int& use_short = m_tinker_mod(rotbnd, use_short);
#endif

} TINKER_NAMESPACE_END

#endif
