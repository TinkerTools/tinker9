#pragma once

#include "macro.h"
#include "sizes.hh"

TINKER_NAMESPACE_BEGIN namespace gkstuf {
using namespace sizes;

extern double& gkc;
extern double (&gkr)[maxtyp];

#ifdef TINKER_MOD_CPP_
extern "C" double TINKER_MOD(gkstuf, gkc);
extern "C" double TINKER_MOD(gkstuf, gkr)[maxtyp];

double& gkc = TINKER_MOD(gkstuf, gkc);
double (&gkr)[maxtyp] = TINKER_MOD(gkstuf, gkr);
#endif
} TINKER_NAMESPACE_END
