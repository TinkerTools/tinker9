#ifndef TINKER_MOD_IELSCF_HH_
#define TINKER_MOD_IELSCF_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace ielscf {
extern int& nfree_aux;
extern double& tautemp_aux;
extern double& kelvin_aux;
extern double*& uaux;
extern double*& upaux;
extern double*& vaux;
extern double*& vpaux;
extern double*& aaux;
extern double*& apaux;
extern int& use_ielscf;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(ielscf, nfree_aux);
extern "C" double m_tinker_mod(ielscf, tautemp_aux);
extern "C" double m_tinker_mod(ielscf, kelvin_aux);
extern "C" double* m_tinker_mod(ielscf, uaux);
extern "C" double* m_tinker_mod(ielscf, upaux);
extern "C" double* m_tinker_mod(ielscf, vaux);
extern "C" double* m_tinker_mod(ielscf, vpaux);
extern "C" double* m_tinker_mod(ielscf, aaux);
extern "C" double* m_tinker_mod(ielscf, apaux);
extern "C" int m_tinker_mod(ielscf, use_ielscf);

int& nfree_aux = m_tinker_mod(ielscf, nfree_aux);
double& tautemp_aux = m_tinker_mod(ielscf, tautemp_aux);
double& kelvin_aux = m_tinker_mod(ielscf, kelvin_aux);
double*& uaux = m_tinker_mod(ielscf, uaux);
double*& upaux = m_tinker_mod(ielscf, upaux);
double*& vaux = m_tinker_mod(ielscf, vaux);
double*& vpaux = m_tinker_mod(ielscf, vpaux);
double*& aaux = m_tinker_mod(ielscf, aaux);
double*& apaux = m_tinker_mod(ielscf, apaux);
int& use_ielscf = m_tinker_mod(ielscf, use_ielscf);
#endif

} TINKER_NAMESPACE_END

#endif
