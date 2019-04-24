#ifndef TINKER_MOD_MUTANT_HH_
#define TINKER_MOD_MUTANT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace mutant {
extern int& nmut;
extern int& vcouple;
extern int*& imut;
extern int*& type0;
extern int*& class0;
extern int*& type1;
extern int*& class1;
extern double& lambda;
extern double& tlambda;
extern double& vlambda;
extern double& elambda;
extern double& scexp;
extern double& scalpha;
extern int*& mut;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(mutant, nmut);
extern "C" int m_tinker_mod(mutant, vcouple);
extern "C" int* m_tinker_mod(mutant, imut);
extern "C" int* m_tinker_mod(mutant, type0);
extern "C" int* m_tinker_mod(mutant, class0);
extern "C" int* m_tinker_mod(mutant, type1);
extern "C" int* m_tinker_mod(mutant, class1);
extern "C" double m_tinker_mod(mutant, lambda);
extern "C" double m_tinker_mod(mutant, tlambda);
extern "C" double m_tinker_mod(mutant, vlambda);
extern "C" double m_tinker_mod(mutant, elambda);
extern "C" double m_tinker_mod(mutant, scexp);
extern "C" double m_tinker_mod(mutant, scalpha);
extern "C" int* m_tinker_mod(mutant, mut);

int& nmut = m_tinker_mod(mutant, nmut);
int& vcouple = m_tinker_mod(mutant, vcouple);
int*& imut = m_tinker_mod(mutant, imut);
int*& type0 = m_tinker_mod(mutant, type0);
int*& class0 = m_tinker_mod(mutant, class0);
int*& type1 = m_tinker_mod(mutant, type1);
int*& class1 = m_tinker_mod(mutant, class1);
double& lambda = m_tinker_mod(mutant, lambda);
double& tlambda = m_tinker_mod(mutant, tlambda);
double& vlambda = m_tinker_mod(mutant, vlambda);
double& elambda = m_tinker_mod(mutant, elambda);
double& scexp = m_tinker_mod(mutant, scexp);
double& scalpha = m_tinker_mod(mutant, scalpha);
int*& mut = m_tinker_mod(mutant, mut);
#endif

} TINKER_NAMESPACE_END

#endif
