#ifndef TINKER_MOD_LINMIN_HH_
#define TINKER_MOD_LINMIN_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace linmin {
extern int& intmax;
extern double& stpmin;
extern double& stpmax;
extern double& cappa;
extern double& slpmax;
extern double& angmax;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(linmin, intmax);
extern "C" double m_tinker_mod(linmin, stpmin);
extern "C" double m_tinker_mod(linmin, stpmax);
extern "C" double m_tinker_mod(linmin, cappa);
extern "C" double m_tinker_mod(linmin, slpmax);
extern "C" double m_tinker_mod(linmin, angmax);

int& intmax = m_tinker_mod(linmin, intmax);
double& stpmin = m_tinker_mod(linmin, stpmin);
double& stpmax = m_tinker_mod(linmin, stpmax);
double& cappa = m_tinker_mod(linmin, cappa);
double& slpmax = m_tinker_mod(linmin, slpmax);
double& angmax = m_tinker_mod(linmin, angmax);
#endif

} TINKER_NAMESPACE_END

#endif
