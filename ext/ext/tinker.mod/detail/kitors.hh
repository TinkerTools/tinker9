#ifndef TINKER_MOD_KITORS_HH_
#define TINKER_MOD_KITORS_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kitors {
const int maxnti = 500;
extern double (&ti1)[maxnti][2];
extern double (&ti2)[maxnti][2];
extern double (&ti3)[maxnti][2];
extern char (&kti)[maxnti][16];

#ifdef TINKER_MOD_CPP_
extern "C" double TINKER_MOD(kitors, ti1)[maxnti][2];
extern "C" double TINKER_MOD(kitors, ti2)[maxnti][2];
extern "C" double TINKER_MOD(kitors, ti3)[maxnti][2];
extern "C" char TINKER_MOD(kitors, kti)[maxnti][16];

double (&ti1)[maxnti][2] = TINKER_MOD(kitors, ti1);
double (&ti2)[maxnti][2] = TINKER_MOD(kitors, ti2);
double (&ti3)[maxnti][2] = TINKER_MOD(kitors, ti3);
char (&kti)[maxnti][16] = TINKER_MOD(kitors, kti);
#endif

} TINKER_NAMESPACE_END

#endif
