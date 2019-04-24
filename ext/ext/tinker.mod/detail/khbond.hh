#ifndef TINKER_MOD_KHBOND_HH_
#define TINKER_MOD_KHBOND_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace khbond {
const int maxnhb = 500;
extern double (&radhb)[maxnhb];
extern double (&epshb)[maxnhb];
extern char (&khb)[maxnhb][8];

#ifdef TINKER_MOD_CPP_
extern "C" double TINKER_MOD(khbond, radhb)[maxnhb];
extern "C" double TINKER_MOD(khbond, epshb)[maxnhb];
extern "C" char TINKER_MOD(khbond, khb)[maxnhb][8];

double (&radhb)[maxnhb] = TINKER_MOD(khbond, radhb);
double (&epshb)[maxnhb] = TINKER_MOD(khbond, epshb);
char (&khb)[maxnhb][8] = TINKER_MOD(khbond, khb);
#endif

} TINKER_NAMESPACE_END

#endif
