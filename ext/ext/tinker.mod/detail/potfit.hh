#ifndef TINKER_MOD_POTFIT_HH_
#define TINKER_MOD_POTFIT_HH_

#include "util/macro.h"
#include "sizes.hh"

TINKER_NAMESPACE_BEGIN namespace potfit {
using namespace sizes;

extern int& nconf;
extern int& namax;
extern int& ngatm;
extern int& nfatm;
extern int (&npgrid)[maxref];
extern int*& ipgrid;
extern double& resp;
extern double (&xdpl0)[maxref];
extern double (&ydpl0)[maxref];
extern double (&zdpl0)[maxref];
extern double (&xxqdp0)[maxref];
extern double (&xyqdp0)[maxref];
extern double (&xzqdp0)[maxref];
extern double (&yyqdp0)[maxref];
extern double (&yzqdp0)[maxref];
extern double (&zzqdp0)[maxref];
extern double*& fit0;
extern double*& fchg;
extern double*& fpol;
extern double*& pgrid;
extern double*& epot;
extern int& use_dpl;
extern int& use_qdp;
extern int& fit_mpl;
extern int& fit_dpl;
extern int& fit_qdp;
extern int*& fitchg;
extern int*& fitpol;
extern int*& gatm;
extern int*& fatm;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(potfit, nconf);
extern "C" int m_tinker_mod(potfit, namax);
extern "C" int m_tinker_mod(potfit, ngatm);
extern "C" int m_tinker_mod(potfit, nfatm);
extern "C" int m_tinker_mod(potfit, npgrid)[maxref];
extern "C" int* m_tinker_mod(potfit, ipgrid);
extern "C" double m_tinker_mod(potfit, resp);
extern "C" double m_tinker_mod(potfit, xdpl0)[maxref];
extern "C" double m_tinker_mod(potfit, ydpl0)[maxref];
extern "C" double m_tinker_mod(potfit, zdpl0)[maxref];
extern "C" double m_tinker_mod(potfit, xxqdp0)[maxref];
extern "C" double m_tinker_mod(potfit, xyqdp0)[maxref];
extern "C" double m_tinker_mod(potfit, xzqdp0)[maxref];
extern "C" double m_tinker_mod(potfit, yyqdp0)[maxref];
extern "C" double m_tinker_mod(potfit, yzqdp0)[maxref];
extern "C" double m_tinker_mod(potfit, zzqdp0)[maxref];
extern "C" double* m_tinker_mod(potfit, fit0);
extern "C" double* m_tinker_mod(potfit, fchg);
extern "C" double* m_tinker_mod(potfit, fpol);
extern "C" double* m_tinker_mod(potfit, pgrid);
extern "C" double* m_tinker_mod(potfit, epot);
extern "C" int m_tinker_mod(potfit, use_dpl);
extern "C" int m_tinker_mod(potfit, use_qdp);
extern "C" int m_tinker_mod(potfit, fit_mpl);
extern "C" int m_tinker_mod(potfit, fit_dpl);
extern "C" int m_tinker_mod(potfit, fit_qdp);
extern "C" int* m_tinker_mod(potfit, fitchg);
extern "C" int* m_tinker_mod(potfit, fitpol);
extern "C" int* m_tinker_mod(potfit, gatm);
extern "C" int* m_tinker_mod(potfit, fatm);

int& nconf = m_tinker_mod(potfit, nconf);
int& namax = m_tinker_mod(potfit, namax);
int& ngatm = m_tinker_mod(potfit, ngatm);
int& nfatm = m_tinker_mod(potfit, nfatm);
int (&npgrid)[maxref] = m_tinker_mod(potfit, npgrid);
int*& ipgrid = m_tinker_mod(potfit, ipgrid);
double& resp = m_tinker_mod(potfit, resp);
double (&xdpl0)[maxref] = m_tinker_mod(potfit, xdpl0);
double (&ydpl0)[maxref] = m_tinker_mod(potfit, ydpl0);
double (&zdpl0)[maxref] = m_tinker_mod(potfit, zdpl0);
double (&xxqdp0)[maxref] = m_tinker_mod(potfit, xxqdp0);
double (&xyqdp0)[maxref] = m_tinker_mod(potfit, xyqdp0);
double (&xzqdp0)[maxref] = m_tinker_mod(potfit, xzqdp0);
double (&yyqdp0)[maxref] = m_tinker_mod(potfit, yyqdp0);
double (&yzqdp0)[maxref] = m_tinker_mod(potfit, yzqdp0);
double (&zzqdp0)[maxref] = m_tinker_mod(potfit, zzqdp0);
double*& fit0 = m_tinker_mod(potfit, fit0);
double*& fchg = m_tinker_mod(potfit, fchg);
double*& fpol = m_tinker_mod(potfit, fpol);
double*& pgrid = m_tinker_mod(potfit, pgrid);
double*& epot = m_tinker_mod(potfit, epot);
int& use_dpl = m_tinker_mod(potfit, use_dpl);
int& use_qdp = m_tinker_mod(potfit, use_qdp);
int& fit_mpl = m_tinker_mod(potfit, fit_mpl);
int& fit_dpl = m_tinker_mod(potfit, fit_dpl);
int& fit_qdp = m_tinker_mod(potfit, fit_qdp);
int*& fitchg = m_tinker_mod(potfit, fitchg);
int*& fitpol = m_tinker_mod(potfit, fitpol);
int*& gatm = m_tinker_mod(potfit, gatm);
int*& fatm = m_tinker_mod(potfit, fatm);
#endif

} TINKER_NAMESPACE_END

#endif
