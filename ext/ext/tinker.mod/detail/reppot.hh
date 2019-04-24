#ifndef TINKER_MOD_REPPOT_HH_
#define TINKER_MOD_REPPOT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace reppot {
extern double& r2scale;
extern double& r3scale;
extern double& r4scale;
extern double& r5scale;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(reppot, r2scale);
extern "C" double m_tinker_mod(reppot, r3scale);
extern "C" double m_tinker_mod(reppot, r4scale);
extern "C" double m_tinker_mod(reppot, r5scale);

double& r2scale = m_tinker_mod(reppot, r2scale);
double& r3scale = m_tinker_mod(reppot, r3scale);
double& r4scale = m_tinker_mod(reppot, r4scale);
double& r5scale = m_tinker_mod(reppot, r5scale);
#endif

} TINKER_NAMESPACE_END

#endif
