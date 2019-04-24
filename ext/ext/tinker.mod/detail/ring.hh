#ifndef TINKER_MOD_RING_HH_
#define TINKER_MOD_RING_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace ring {
extern int& nring3;
extern int& nring4;
extern int& nring5;
extern int& nring6;
extern int*& iring3;
extern int*& iring4;
extern int*& iring5;
extern int*& iring6;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(ring, nring3);
extern "C" int m_tinker_mod(ring, nring4);
extern "C" int m_tinker_mod(ring, nring5);
extern "C" int m_tinker_mod(ring, nring6);
extern "C" int* m_tinker_mod(ring, iring3);
extern "C" int* m_tinker_mod(ring, iring4);
extern "C" int* m_tinker_mod(ring, iring5);
extern "C" int* m_tinker_mod(ring, iring6);

int& nring3 = m_tinker_mod(ring, nring3);
int& nring4 = m_tinker_mod(ring, nring4);
int& nring5 = m_tinker_mod(ring, nring5);
int& nring6 = m_tinker_mod(ring, nring6);
int*& iring3 = m_tinker_mod(ring, iring3);
int*& iring4 = m_tinker_mod(ring, iring4);
int*& iring5 = m_tinker_mod(ring, iring5);
int*& iring6 = m_tinker_mod(ring, iring6);
#endif

} TINKER_NAMESPACE_END

#endif
