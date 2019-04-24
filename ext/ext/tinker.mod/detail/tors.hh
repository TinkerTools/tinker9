#ifndef TINKER_MOD_TORS_HH_
#define TINKER_MOD_TORS_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace tors {
extern int& ntors;
extern int*& itors;
extern double*& tors1;
extern double*& tors2;
extern double*& tors3;
extern double*& tors4;
extern double*& tors5;
extern double*& tors6;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(tors, ntors);
extern "C" int* m_tinker_mod(tors, itors);
extern "C" double* m_tinker_mod(tors, tors1);
extern "C" double* m_tinker_mod(tors, tors2);
extern "C" double* m_tinker_mod(tors, tors3);
extern "C" double* m_tinker_mod(tors, tors4);
extern "C" double* m_tinker_mod(tors, tors5);
extern "C" double* m_tinker_mod(tors, tors6);

int& ntors = m_tinker_mod(tors, ntors);
int*& itors = m_tinker_mod(tors, itors);
double*& tors1 = m_tinker_mod(tors, tors1);
double*& tors2 = m_tinker_mod(tors, tors2);
double*& tors3 = m_tinker_mod(tors, tors3);
double*& tors4 = m_tinker_mod(tors, tors4);
double*& tors5 = m_tinker_mod(tors, tors5);
double*& tors6 = m_tinker_mod(tors, tors6);
#endif

} TINKER_NAMESPACE_END

#endif
