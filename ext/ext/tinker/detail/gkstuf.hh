#pragma once

#include "macro.h"
#include "sizes.hh"

namespace tinker { namespace gkstuf {
using namespace sizes;

extern double& gkc;
extern double (&gkr)[maxtyp];

#ifdef TINKER_FORTRAN_MODULE_CPP
extern "C" double TINKER_MOD(gkstuf, gkc);
extern "C" double TINKER_MOD(gkstuf, gkr)[maxtyp];

double& gkc = TINKER_MOD(gkstuf, gkc);
double (&gkr)[maxtyp] = TINKER_MOD(gkstuf, gkr);
#endif
} }
