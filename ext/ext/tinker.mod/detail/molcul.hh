#ifndef TINKER_MOD_MOLCUL_HH_
#define TINKER_MOD_MOLCUL_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace molcul {
extern int& nmol;
extern int*& imol;
extern int*& kmol;
extern int*& molcule;
extern double& totmass;
extern double*& molmass;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(molcul, nmol);
extern "C" int* m_tinker_mod(molcul, imol);
extern "C" int* m_tinker_mod(molcul, kmol);
extern "C" int* m_tinker_mod(molcul, molcule);
extern "C" double m_tinker_mod(molcul, totmass);
extern "C" double* m_tinker_mod(molcul, molmass);

int& nmol = m_tinker_mod(molcul, nmol);
int*& imol = m_tinker_mod(molcul, imol);
int*& kmol = m_tinker_mod(molcul, kmol);
int*& molcule = m_tinker_mod(molcul, molcule);
double& totmass = m_tinker_mod(molcul, totmass);
double*& molmass = m_tinker_mod(molcul, molmass);
#endif

} TINKER_NAMESPACE_END

#endif
