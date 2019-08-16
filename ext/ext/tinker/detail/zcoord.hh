#ifndef TINKER_MOD_ZCOORD_HH_
#define TINKER_MOD_ZCOORD_HH_

#include "macro.h"
#include "sizes.hh"

TINKER_NAMESPACE_BEGIN namespace zcoord {
using namespace sizes;

extern int (&iz)[maxatm][4];
extern double (&zbond)[maxatm];
extern double (&zang)[maxatm];
extern double (&ztors)[maxatm];

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(zcoord, iz)[maxatm][4];
extern "C" double TINKER_MOD(zcoord, zbond)[maxatm];
extern "C" double TINKER_MOD(zcoord, zang)[maxatm];
extern "C" double TINKER_MOD(zcoord, ztors)[maxatm];

int (&iz)[maxatm][4] = TINKER_MOD(zcoord, iz);
double (&zbond)[maxatm] = TINKER_MOD(zcoord, zbond);
double (&zang)[maxatm] = TINKER_MOD(zcoord, zang);
double (&ztors)[maxatm] = TINKER_MOD(zcoord, ztors);
#endif
} TINKER_NAMESPACE_END

#endif
