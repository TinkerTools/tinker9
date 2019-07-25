#ifndef TINKER_MOD_BOUND_HH_
#define TINKER_MOD_BOUND_HH_

#include "util_macro.h"

TINKER_NAMESPACE_BEGIN namespace bound {
extern double& polycut;
extern double& polycut2;
extern int& use_bounds;
extern int& use_replica;
extern int& use_polymer;

#ifdef TINKER_MOD_CPP_
extern "C" double TINKER_MOD(bound, polycut);
extern "C" double TINKER_MOD(bound, polycut2);
extern "C" int TINKER_MOD(bound, use_bounds);
extern "C" int TINKER_MOD(bound, use_replica);
extern "C" int TINKER_MOD(bound, use_polymer);

double& polycut = TINKER_MOD(bound, polycut);
double& polycut2 = TINKER_MOD(bound, polycut2);
int& use_bounds = TINKER_MOD(bound, use_bounds);
int& use_replica = TINKER_MOD(bound, use_replica);
int& use_polymer = TINKER_MOD(bound, use_polymer);
#endif
} TINKER_NAMESPACE_END

#endif
