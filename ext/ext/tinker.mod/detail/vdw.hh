#ifndef TINKER_MOD_VDW_HH_
#define TINKER_MOD_VDW_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace vdw {
extern int& nvdw;
extern int*& ivdw;
extern int*& jvdw;
extern int*& ired;
extern double*& kred;
extern double*& xred;
extern double*& yred;
extern double*& zred;
extern double*& radmin;
extern double*& epsilon;
extern double*& radmin4;
extern double*& epsilon4;
extern double*& radhbnd;
extern double*& epshbnd;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(vdw, nvdw);
extern "C" int* m_tinker_mod(vdw, ivdw);
extern "C" int* m_tinker_mod(vdw, jvdw);
extern "C" int* m_tinker_mod(vdw, ired);
extern "C" double* m_tinker_mod(vdw, kred);
extern "C" double* m_tinker_mod(vdw, xred);
extern "C" double* m_tinker_mod(vdw, yred);
extern "C" double* m_tinker_mod(vdw, zred);
extern "C" double* m_tinker_mod(vdw, radmin);
extern "C" double* m_tinker_mod(vdw, epsilon);
extern "C" double* m_tinker_mod(vdw, radmin4);
extern "C" double* m_tinker_mod(vdw, epsilon4);
extern "C" double* m_tinker_mod(vdw, radhbnd);
extern "C" double* m_tinker_mod(vdw, epshbnd);

int& nvdw = m_tinker_mod(vdw, nvdw);
int*& ivdw = m_tinker_mod(vdw, ivdw);
int*& jvdw = m_tinker_mod(vdw, jvdw);
int*& ired = m_tinker_mod(vdw, ired);
double*& kred = m_tinker_mod(vdw, kred);
double*& xred = m_tinker_mod(vdw, xred);
double*& yred = m_tinker_mod(vdw, yred);
double*& zred = m_tinker_mod(vdw, zred);
double*& radmin = m_tinker_mod(vdw, radmin);
double*& epsilon = m_tinker_mod(vdw, epsilon);
double*& radmin4 = m_tinker_mod(vdw, radmin4);
double*& epsilon4 = m_tinker_mod(vdw, epsilon4);
double*& radhbnd = m_tinker_mod(vdw, radhbnd);
double*& epshbnd = m_tinker_mod(vdw, epshbnd);
#endif

} TINKER_NAMESPACE_END

#endif
