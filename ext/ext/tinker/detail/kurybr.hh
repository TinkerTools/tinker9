#ifndef TINKER_MOD_KURYBR_HH_
#define TINKER_MOD_KURYBR_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kurybr {
const int maxnu = 2000;
extern double (&ucon)[maxnu];
extern double (&dst13)[maxnu];
extern char (&ku)[maxnu][12];

#ifdef TINKER_MOD_CPP_
extern "C" double TINKER_MOD(kurybr, ucon)[maxnu];
extern "C" double TINKER_MOD(kurybr, dst13)[maxnu];
extern "C" char TINKER_MOD(kurybr, ku)[maxnu][12];

double (&ucon)[maxnu] = TINKER_MOD(kurybr, ucon);
double (&dst13)[maxnu] = TINKER_MOD(kurybr, dst13);
char (&ku)[maxnu][12] = TINKER_MOD(kurybr, ku);
#endif
} TINKER_NAMESPACE_END

#endif
