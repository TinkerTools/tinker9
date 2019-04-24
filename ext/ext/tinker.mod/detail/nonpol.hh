#ifndef TINKER_MOD_NONPOL_HH_
#define TINKER_MOD_NONPOL_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace nonpol {
const double epso = 0.1100e0;
const double epsh = 0.0135e0;
const double rmino = 1.7025e0;
const double rminh = 1.3275e0;
const double awater = 0.033428e0;
const double slevy = 1.0e0;
extern double& solvprs;
extern double& surften;
extern double& spcut;
extern double& spoff;
extern double& stcut;
extern double& stoff;
extern double*& rcav;
extern double*& rdisp;
extern double*& cdisp;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(nonpol, solvprs);
extern "C" double m_tinker_mod(nonpol, surften);
extern "C" double m_tinker_mod(nonpol, spcut);
extern "C" double m_tinker_mod(nonpol, spoff);
extern "C" double m_tinker_mod(nonpol, stcut);
extern "C" double m_tinker_mod(nonpol, stoff);
extern "C" double* m_tinker_mod(nonpol, rcav);
extern "C" double* m_tinker_mod(nonpol, rdisp);
extern "C" double* m_tinker_mod(nonpol, cdisp);

double& solvprs = m_tinker_mod(nonpol, solvprs);
double& surften = m_tinker_mod(nonpol, surften);
double& spcut = m_tinker_mod(nonpol, spcut);
double& spoff = m_tinker_mod(nonpol, spoff);
double& stcut = m_tinker_mod(nonpol, stcut);
double& stoff = m_tinker_mod(nonpol, stoff);
double*& rcav = m_tinker_mod(nonpol, rcav);
double*& rdisp = m_tinker_mod(nonpol, rdisp);
double*& cdisp = m_tinker_mod(nonpol, cdisp);
#endif

} TINKER_NAMESPACE_END

#endif
