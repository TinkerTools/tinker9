#ifndef TINKER_MOD_QMSTUF_HH_
#define TINKER_MOD_QMSTUF_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace qmstuf {
extern int& ngatom;
extern double& egau;
extern double*& gx;
extern double*& gy;
extern double*& gz;
extern double*& gfreq;
extern double*& gforce;
extern double*& gh;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(qmstuf, ngatom);
extern "C" double m_tinker_mod(qmstuf, egau);
extern "C" double* m_tinker_mod(qmstuf, gx);
extern "C" double* m_tinker_mod(qmstuf, gy);
extern "C" double* m_tinker_mod(qmstuf, gz);
extern "C" double* m_tinker_mod(qmstuf, gfreq);
extern "C" double* m_tinker_mod(qmstuf, gforce);
extern "C" double* m_tinker_mod(qmstuf, gh);

int& ngatom = m_tinker_mod(qmstuf, ngatom);
double& egau = m_tinker_mod(qmstuf, egau);
double*& gx = m_tinker_mod(qmstuf, gx);
double*& gy = m_tinker_mod(qmstuf, gy);
double*& gz = m_tinker_mod(qmstuf, gz);
double*& gfreq = m_tinker_mod(qmstuf, gfreq);
double*& gforce = m_tinker_mod(qmstuf, gforce);
double*& gh = m_tinker_mod(qmstuf, gh);
#endif

} TINKER_NAMESPACE_END

#endif
