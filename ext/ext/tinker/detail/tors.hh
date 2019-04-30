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
extern "C" int TINKER_MOD(tors, ntors);
extern "C" int* TINKER_MOD(tors, itors);
extern "C" double* TINKER_MOD(tors, tors1);
extern "C" double* TINKER_MOD(tors, tors2);
extern "C" double* TINKER_MOD(tors, tors3);
extern "C" double* TINKER_MOD(tors, tors4);
extern "C" double* TINKER_MOD(tors, tors5);
extern "C" double* TINKER_MOD(tors, tors6);

int& ntors = TINKER_MOD(tors, ntors);
int*& itors = TINKER_MOD(tors, itors);
double*& tors1 = TINKER_MOD(tors, tors1);
double*& tors2 = TINKER_MOD(tors, tors2);
double*& tors3 = TINKER_MOD(tors, tors3);
double*& tors4 = TINKER_MOD(tors, tors4);
double*& tors5 = TINKER_MOD(tors, tors5);
double*& tors6 = TINKER_MOD(tors, tors6);
#endif

} TINKER_NAMESPACE_END

#endif
