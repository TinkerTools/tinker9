#ifndef TINKER_MOD_BOUND_HH_
#define TINKER_MOD_BOUND_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace bound {
extern double& polycut;
extern double& polycut2;
extern int& use_bounds;
extern int& use_replica;
extern int& use_polymer;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(bound, polycut);
extern "C" double m_tinker_mod(bound, polycut2);
extern "C" int m_tinker_mod(bound, use_bounds);
extern "C" int m_tinker_mod(bound, use_replica);
extern "C" int m_tinker_mod(bound, use_polymer);

double& polycut = m_tinker_mod(bound, polycut);
double& polycut2 = m_tinker_mod(bound, polycut2);
int& use_bounds = m_tinker_mod(bound, use_bounds);
int& use_replica = m_tinker_mod(bound, use_replica);
int& use_polymer = m_tinker_mod(bound, use_polymer);
#endif

} TINKER_NAMESPACE_END

#endif
