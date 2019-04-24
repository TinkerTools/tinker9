#ifndef TINKER_MOD_KOPDST_HH_
#define TINKER_MOD_KOPDST_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kopdst {
const int maxnopd = 500;
extern double (&opds)[maxnopd];
extern char (&kopd)[maxnopd][16];

#ifdef TINKER_MOD_CPP_
extern "C" double TINKER_MOD(kopdst, opds)[maxnopd];
extern "C" char TINKER_MOD(kopdst, kopd)[maxnopd][16];

double (&opds)[maxnopd] = TINKER_MOD(kopdst, opds);
char (&kopd)[maxnopd][16] = TINKER_MOD(kopdst, kopd);
#endif

} TINKER_NAMESPACE_END

#endif
