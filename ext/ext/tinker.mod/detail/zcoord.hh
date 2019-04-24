#ifndef TINKER_MOD_ZCOORD_HH_
#define TINKER_MOD_ZCOORD_HH_

#include "util/macro.h"
#include "sizes.hh"

TINKER_NAMESPACE_BEGIN namespace zcoord {
using namespace sizes;

extern int (&iz)[maxatm][4];
extern double (&zbond)[maxatm];
extern double (&zang)[maxatm];
extern double (&ztors)[maxatm];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(zcoord, iz)[maxatm][4];
extern "C" double m_tinker_mod(zcoord, zbond)[maxatm];
extern "C" double m_tinker_mod(zcoord, zang)[maxatm];
extern "C" double m_tinker_mod(zcoord, ztors)[maxatm];

int (&iz)[maxatm][4] = m_tinker_mod(zcoord, iz);
double (&zbond)[maxatm] = m_tinker_mod(zcoord, zbond);
double (&zang)[maxatm] = m_tinker_mod(zcoord, zang);
double (&ztors)[maxatm] = m_tinker_mod(zcoord, ztors);
#endif

} TINKER_NAMESPACE_END

#endif
