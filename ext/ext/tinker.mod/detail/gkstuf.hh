#ifndef TINKER_MOD_GKSTUF_HH_
#define TINKER_MOD_GKSTUF_HH_

#include "util/macro.h"
#include "sizes.hh"

TINKER_NAMESPACE_BEGIN namespace gkstuf {
using namespace sizes;

extern double& gkc;
extern double (&gkr)[maxtyp];

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(gkstuf, gkc);
extern "C" double m_tinker_mod(gkstuf, gkr)[maxtyp];

double& gkc = m_tinker_mod(gkstuf, gkc);
double (&gkr)[maxtyp] = m_tinker_mod(gkstuf, gkr);
#endif

} TINKER_NAMESPACE_END

#endif
