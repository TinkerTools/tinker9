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
extern "C" int TINKER_MOD(potfit, nconf);
extern "C" int TINKER_MOD(potfit, namax);
extern "C" int TINKER_MOD(potfit, ngatm);
extern "C" int TINKER_MOD(potfit, nfatm);
extern "C" int TINKER_MOD(potfit, npgrid)[maxref];
extern "C" int* TINKER_MOD(potfit, ipgrid);
extern "C" double TINKER_MOD(potfit, resp);
extern "C" double TINKER_MOD(potfit, xdpl0)[maxref];
extern "C" double TINKER_MOD(potfit, ydpl0)[maxref];
extern "C" double TINKER_MOD(potfit, zdpl0)[maxref];
extern "C" double TINKER_MOD(potfit, xxqdp0)[maxref];
extern "C" double TINKER_MOD(potfit, xyqdp0)[maxref];
extern "C" double TINKER_MOD(potfit, xzqdp0)[maxref];
extern "C" double TINKER_MOD(potfit, yyqdp0)[maxref];
extern "C" double TINKER_MOD(potfit, yzqdp0)[maxref];
extern "C" double TINKER_MOD(potfit, zzqdp0)[maxref];
extern "C" double* TINKER_MOD(potfit, fit0);
extern "C" double* TINKER_MOD(potfit, fchg);
extern "C" double* TINKER_MOD(potfit, fpol);
extern "C" double* TINKER_MOD(potfit, pgrid);
extern "C" double* TINKER_MOD(potfit, epot);
extern "C" int TINKER_MOD(potfit, use_dpl);
extern "C" int TINKER_MOD(potfit, use_qdp);
extern "C" int TINKER_MOD(potfit, fit_mpl);
extern "C" int TINKER_MOD(potfit, fit_dpl);
extern "C" int TINKER_MOD(potfit, fit_qdp);
extern "C" int* TINKER_MOD(potfit, fitchg);
extern "C" int* TINKER_MOD(potfit, fitpol);
extern "C" int* TINKER_MOD(potfit, gatm);
extern "C" int* TINKER_MOD(potfit, fatm);

int& nconf = TINKER_MOD(potfit, nconf);
int& namax = TINKER_MOD(potfit, namax);
int& ngatm = TINKER_MOD(potfit, ngatm);
int& nfatm = TINKER_MOD(potfit, nfatm);
int (&npgrid)[maxref] = TINKER_MOD(potfit, npgrid);
int*& ipgrid = TINKER_MOD(potfit, ipgrid);
double& resp = TINKER_MOD(potfit, resp);
double (&xdpl0)[maxref] = TINKER_MOD(potfit, xdpl0);
double (&ydpl0)[maxref] = TINKER_MOD(potfit, ydpl0);
double (&zdpl0)[maxref] = TINKER_MOD(potfit, zdpl0);
double (&xxqdp0)[maxref] = TINKER_MOD(potfit, xxqdp0);
double (&xyqdp0)[maxref] = TINKER_MOD(potfit, xyqdp0);
double (&xzqdp0)[maxref] = TINKER_MOD(potfit, xzqdp0);
double (&yyqdp0)[maxref] = TINKER_MOD(potfit, yyqdp0);
double (&yzqdp0)[maxref] = TINKER_MOD(potfit, yzqdp0);
double (&zzqdp0)[maxref] = TINKER_MOD(potfit, zzqdp0);
double*& fit0 = TINKER_MOD(potfit, fit0);
double*& fchg = TINKER_MOD(potfit, fchg);
double*& fpol = TINKER_MOD(potfit, fpol);
double*& pgrid = TINKER_MOD(potfit, pgrid);
double*& epot = TINKER_MOD(potfit, epot);
int& use_dpl = TINKER_MOD(potfit, use_dpl);
int& use_qdp = TINKER_MOD(potfit, use_qdp);
int& fit_mpl = TINKER_MOD(potfit, fit_mpl);
int& fit_dpl = TINKER_MOD(potfit, fit_dpl);
int& fit_qdp = TINKER_MOD(potfit, fit_qdp);
int*& fitchg = TINKER_MOD(potfit, fitchg);
int*& fitpol = TINKER_MOD(potfit, fitpol);
int*& gatm = TINKER_MOD(potfit, gatm);
int*& fatm = TINKER_MOD(potfit, fatm);
#endif
} TINKER_NAMESPACE_END

#endif
