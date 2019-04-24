#ifndef TINKER_MOD_MPLPOT_HH_
#define TINKER_MOD_MPLPOT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace mplpot {
extern double& m2scale;
extern double& m3scale;
extern double& m4scale;
extern double& m5scale;
extern int& use_chgpen;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(mplpot, m2scale);
extern "C" double m_tinker_mod(mplpot, m3scale);
extern "C" double m_tinker_mod(mplpot, m4scale);
extern "C" double m_tinker_mod(mplpot, m5scale);
extern "C" int m_tinker_mod(mplpot, use_chgpen);

double& m2scale = m_tinker_mod(mplpot, m2scale);
double& m3scale = m_tinker_mod(mplpot, m3scale);
double& m4scale = m_tinker_mod(mplpot, m4scale);
double& m5scale = m_tinker_mod(mplpot, m5scale);
int& use_chgpen = m_tinker_mod(mplpot, use_chgpen);
#endif

} TINKER_NAMESPACE_END

#endif
