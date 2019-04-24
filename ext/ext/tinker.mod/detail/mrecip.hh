#ifndef TINKER_MOD_MRECIP_HH_
#define TINKER_MOD_MRECIP_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace mrecip {
extern double& vmxx;
extern double& vmyy;
extern double& vmzz;
extern double& vmxy;
extern double& vmxz;
extern double& vmyz;
extern double*& cmp;
extern double*& fmp;
extern double*& cphi;
extern double*& fphi;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(mrecip, vmxx);
extern "C" double m_tinker_mod(mrecip, vmyy);
extern "C" double m_tinker_mod(mrecip, vmzz);
extern "C" double m_tinker_mod(mrecip, vmxy);
extern "C" double m_tinker_mod(mrecip, vmxz);
extern "C" double m_tinker_mod(mrecip, vmyz);
extern "C" double* m_tinker_mod(mrecip, cmp);
extern "C" double* m_tinker_mod(mrecip, fmp);
extern "C" double* m_tinker_mod(mrecip, cphi);
extern "C" double* m_tinker_mod(mrecip, fphi);

double& vmxx = m_tinker_mod(mrecip, vmxx);
double& vmyy = m_tinker_mod(mrecip, vmyy);
double& vmzz = m_tinker_mod(mrecip, vmzz);
double& vmxy = m_tinker_mod(mrecip, vmxy);
double& vmxz = m_tinker_mod(mrecip, vmxz);
double& vmyz = m_tinker_mod(mrecip, vmyz);
double*& cmp = m_tinker_mod(mrecip, cmp);
double*& fmp = m_tinker_mod(mrecip, fmp);
double*& cphi = m_tinker_mod(mrecip, cphi);
double*& fphi = m_tinker_mod(mrecip, fphi);
#endif

} TINKER_NAMESPACE_END

#endif
