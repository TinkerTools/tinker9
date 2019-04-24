#ifndef TINKER_MOD_POLOPT_HH_
#define TINKER_MOD_POLOPT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace polopt {
const int maxopt = 6;
extern int& optorder;
extern int& optlevel;
extern double*& copt;
extern double*& copm;
extern double*& uopt;
extern double*& uoptp;
extern double*& uopts;
extern double*& uoptps;
extern double*& fopt;
extern double*& foptp;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(polopt, optorder);
extern "C" int m_tinker_mod(polopt, optlevel);
extern "C" double* m_tinker_mod(polopt, copt);
extern "C" double* m_tinker_mod(polopt, copm);
extern "C" double* m_tinker_mod(polopt, uopt);
extern "C" double* m_tinker_mod(polopt, uoptp);
extern "C" double* m_tinker_mod(polopt, uopts);
extern "C" double* m_tinker_mod(polopt, uoptps);
extern "C" double* m_tinker_mod(polopt, fopt);
extern "C" double* m_tinker_mod(polopt, foptp);

int& optorder = m_tinker_mod(polopt, optorder);
int& optlevel = m_tinker_mod(polopt, optlevel);
double*& copt = m_tinker_mod(polopt, copt);
double*& copm = m_tinker_mod(polopt, copm);
double*& uopt = m_tinker_mod(polopt, uopt);
double*& uoptp = m_tinker_mod(polopt, uoptp);
double*& uopts = m_tinker_mod(polopt, uopts);
double*& uoptps = m_tinker_mod(polopt, uoptps);
double*& fopt = m_tinker_mod(polopt, fopt);
double*& foptp = m_tinker_mod(polopt, foptp);
#endif

} TINKER_NAMESPACE_END

#endif
