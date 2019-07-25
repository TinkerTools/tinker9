#ifndef TINKER_MOD_REPPOT_HH_
#define TINKER_MOD_REPPOT_HH_

#include "util_macro.h"

TINKER_NAMESPACE_BEGIN namespace reppot {
extern double& r2scale;
extern double& r3scale;
extern double& r4scale;
extern double& r5scale;

#ifdef TINKER_MOD_CPP_
extern "C" double TINKER_MOD(reppot, r2scale);
extern "C" double TINKER_MOD(reppot, r3scale);
extern "C" double TINKER_MOD(reppot, r4scale);
extern "C" double TINKER_MOD(reppot, r5scale);

double& r2scale = TINKER_MOD(reppot, r2scale);
double& r3scale = TINKER_MOD(reppot, r3scale);
double& r4scale = TINKER_MOD(reppot, r4scale);
double& r5scale = TINKER_MOD(reppot, r5scale);
#endif
} TINKER_NAMESPACE_END

#endif
