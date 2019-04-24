#ifndef TINKER_MOD_RGDDYN_HH_
#define TINKER_MOD_RGDDYN_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace rgddyn {
extern double*& xcmo;
extern double*& ycmo;
extern double*& zcmo;
extern double*& vcm;
extern double*& wcm;
extern double*& lm;
extern double*& vc;
extern double*& wc;
extern int*& linear;

#ifdef TINKER_MOD_CPP_
extern "C" double* m_tinker_mod(rgddyn, xcmo);
extern "C" double* m_tinker_mod(rgddyn, ycmo);
extern "C" double* m_tinker_mod(rgddyn, zcmo);
extern "C" double* m_tinker_mod(rgddyn, vcm);
extern "C" double* m_tinker_mod(rgddyn, wcm);
extern "C" double* m_tinker_mod(rgddyn, lm);
extern "C" double* m_tinker_mod(rgddyn, vc);
extern "C" double* m_tinker_mod(rgddyn, wc);
extern "C" int* m_tinker_mod(rgddyn, linear);

double*& xcmo = m_tinker_mod(rgddyn, xcmo);
double*& ycmo = m_tinker_mod(rgddyn, ycmo);
double*& zcmo = m_tinker_mod(rgddyn, zcmo);
double*& vcm = m_tinker_mod(rgddyn, vcm);
double*& wcm = m_tinker_mod(rgddyn, wcm);
double*& lm = m_tinker_mod(rgddyn, lm);
double*& vc = m_tinker_mod(rgddyn, vc);
double*& wc = m_tinker_mod(rgddyn, wc);
int*& linear = m_tinker_mod(rgddyn, linear);
#endif

} TINKER_NAMESPACE_END

#endif
