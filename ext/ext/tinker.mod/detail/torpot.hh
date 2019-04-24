#ifndef TINKER_MOD_TORPOT_HH_
#define TINKER_MOD_TORPOT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace torpot {
extern double& idihunit;
extern double& itorunit;
extern double& torsunit;
extern double& ptorunit;
extern double& storunit;
extern double& atorunit;
extern double& ttorunit;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(torpot, idihunit);
extern "C" double m_tinker_mod(torpot, itorunit);
extern "C" double m_tinker_mod(torpot, torsunit);
extern "C" double m_tinker_mod(torpot, ptorunit);
extern "C" double m_tinker_mod(torpot, storunit);
extern "C" double m_tinker_mod(torpot, atorunit);
extern "C" double m_tinker_mod(torpot, ttorunit);

double& idihunit = m_tinker_mod(torpot, idihunit);
double& itorunit = m_tinker_mod(torpot, itorunit);
double& torsunit = m_tinker_mod(torpot, torsunit);
double& ptorunit = m_tinker_mod(torpot, ptorunit);
double& storunit = m_tinker_mod(torpot, storunit);
double& atorunit = m_tinker_mod(torpot, atorunit);
double& ttorunit = m_tinker_mod(torpot, ttorunit);
#endif

} TINKER_NAMESPACE_END

#endif
