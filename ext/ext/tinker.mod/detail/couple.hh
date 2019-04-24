#ifndef TINKER_MOD_COUPLE_HH_
#define TINKER_MOD_COUPLE_HH_

#include "util/macro.h"
#include "sizes.hh"

TINKER_NAMESPACE_BEGIN namespace couple {
using namespace sizes;

extern int (&n12)[maxatm];
extern int*& n13;
extern int*& n14;
extern int*& n15;
extern int (&i12)[maxatm][maxval];
extern int*& i13;
extern int*& i14;
extern int*& i15;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(couple, n12)[maxatm];
extern "C" int* m_tinker_mod(couple, n13);
extern "C" int* m_tinker_mod(couple, n14);
extern "C" int* m_tinker_mod(couple, n15);
extern "C" int m_tinker_mod(couple, i12)[maxatm][maxval];
extern "C" int* m_tinker_mod(couple, i13);
extern "C" int* m_tinker_mod(couple, i14);
extern "C" int* m_tinker_mod(couple, i15);

int (&n12)[maxatm] = m_tinker_mod(couple, n12);
int*& n13 = m_tinker_mod(couple, n13);
int*& n14 = m_tinker_mod(couple, n14);
int*& n15 = m_tinker_mod(couple, n15);
int (&i12)[maxatm][maxval] = m_tinker_mod(couple, i12);
int*& i13 = m_tinker_mod(couple, i13);
int*& i14 = m_tinker_mod(couple, i14);
int*& i15 = m_tinker_mod(couple, i15);
#endif

} TINKER_NAMESPACE_END

#endif
