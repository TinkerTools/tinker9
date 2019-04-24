#ifndef TINKER_MOD_PATHS_HH_
#define TINKER_MOD_PATHS_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace paths {
extern double& pnorm;
extern double (&acoeff)[7][7];
extern double*& pc0;
extern double*& pc1;
extern double*& pvect;
extern double*& pstep;
extern double*& pzet;
extern double*& gc;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(paths, pnorm);
extern "C" double m_tinker_mod(paths, acoeff)[7][7];
extern "C" double* m_tinker_mod(paths, pc0);
extern "C" double* m_tinker_mod(paths, pc1);
extern "C" double* m_tinker_mod(paths, pvect);
extern "C" double* m_tinker_mod(paths, pstep);
extern "C" double* m_tinker_mod(paths, pzet);
extern "C" double* m_tinker_mod(paths, gc);

double& pnorm = m_tinker_mod(paths, pnorm);
double (&acoeff)[7][7] = m_tinker_mod(paths, acoeff);
double*& pc0 = m_tinker_mod(paths, pc0);
double*& pc1 = m_tinker_mod(paths, pc1);
double*& pvect = m_tinker_mod(paths, pvect);
double*& pstep = m_tinker_mod(paths, pstep);
double*& pzet = m_tinker_mod(paths, pzet);
double*& gc = m_tinker_mod(paths, gc);
#endif

} TINKER_NAMESPACE_END

#endif
