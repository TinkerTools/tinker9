#ifndef TINKER_MOD_ATMLST_HH_
#define TINKER_MOD_ATMLST_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace atmlst {
extern int*& bndlist;
extern int*& anglist;

#ifdef TINKER_MOD_CPP_
extern "C" int* m_tinker_mod(atmlst, bndlist);
extern "C" int* m_tinker_mod(atmlst, anglist);

int*& bndlist = m_tinker_mod(atmlst, bndlist);
int*& anglist = m_tinker_mod(atmlst, anglist);
#endif

} TINKER_NAMESPACE_END

#endif
