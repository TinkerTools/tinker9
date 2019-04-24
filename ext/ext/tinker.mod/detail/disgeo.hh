#ifndef TINKER_MOD_DISGEO_HH_
#define TINKER_MOD_DISGEO_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace disgeo {
extern double& vdwmax;
extern double& compact;
extern double& pathmax;
extern double*& dbnd;
extern double*& georad;
extern int& use_invert;
extern int& use_anneal;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(disgeo, vdwmax);
extern "C" double m_tinker_mod(disgeo, compact);
extern "C" double m_tinker_mod(disgeo, pathmax);
extern "C" double* m_tinker_mod(disgeo, dbnd);
extern "C" double* m_tinker_mod(disgeo, georad);
extern "C" int m_tinker_mod(disgeo, use_invert);
extern "C" int m_tinker_mod(disgeo, use_anneal);

double& vdwmax = m_tinker_mod(disgeo, vdwmax);
double& compact = m_tinker_mod(disgeo, compact);
double& pathmax = m_tinker_mod(disgeo, pathmax);
double*& dbnd = m_tinker_mod(disgeo, dbnd);
double*& georad = m_tinker_mod(disgeo, georad);
int& use_invert = m_tinker_mod(disgeo, use_invert);
int& use_anneal = m_tinker_mod(disgeo, use_anneal);
#endif

} TINKER_NAMESPACE_END

#endif
