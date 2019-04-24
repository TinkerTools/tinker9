#ifndef TINKER_MOD_DISP_HH_
#define TINKER_MOD_DISP_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace disp {
extern int& ndisp;
extern int*& idisp;
extern double& csixpr;
extern double*& csix;
extern double*& adisp;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(disp, ndisp);
extern "C" int* m_tinker_mod(disp, idisp);
extern "C" double m_tinker_mod(disp, csixpr);
extern "C" double* m_tinker_mod(disp, csix);
extern "C" double* m_tinker_mod(disp, adisp);

int& ndisp = m_tinker_mod(disp, ndisp);
int*& idisp = m_tinker_mod(disp, idisp);
double& csixpr = m_tinker_mod(disp, csixpr);
double*& csix = m_tinker_mod(disp, csix);
double*& adisp = m_tinker_mod(disp, adisp);
#endif

} TINKER_NAMESPACE_END

#endif
