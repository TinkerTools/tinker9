#ifndef TINKER_MOD_MOMENT_HH_
#define TINKER_MOD_MOMENT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace moment {
extern double& netchg;
extern double& netdpl;
extern double (&netqdp)[3];
extern double& xdpl;
extern double& ydpl;
extern double& zdpl;
extern double& xxqdp;
extern double& xyqdp;
extern double& xzqdp;
extern double& yxqdp;
extern double& yyqdp;
extern double& yzqdp;
extern double& zxqdp;
extern double& zyqdp;
extern double& zzqdp;

#ifdef TINKER_MOD_CPP_
extern "C" double TINKER_MOD(moment, netchg);
extern "C" double TINKER_MOD(moment, netdpl);
extern "C" double TINKER_MOD(moment, netqdp)[3];
extern "C" double TINKER_MOD(moment, xdpl);
extern "C" double TINKER_MOD(moment, ydpl);
extern "C" double TINKER_MOD(moment, zdpl);
extern "C" double TINKER_MOD(moment, xxqdp);
extern "C" double TINKER_MOD(moment, xyqdp);
extern "C" double TINKER_MOD(moment, xzqdp);
extern "C" double TINKER_MOD(moment, yxqdp);
extern "C" double TINKER_MOD(moment, yyqdp);
extern "C" double TINKER_MOD(moment, yzqdp);
extern "C" double TINKER_MOD(moment, zxqdp);
extern "C" double TINKER_MOD(moment, zyqdp);
extern "C" double TINKER_MOD(moment, zzqdp);

double& netchg = TINKER_MOD(moment, netchg);
double& netdpl = TINKER_MOD(moment, netdpl);
double (&netqdp)[3] = TINKER_MOD(moment, netqdp);
double& xdpl = TINKER_MOD(moment, xdpl);
double& ydpl = TINKER_MOD(moment, ydpl);
double& zdpl = TINKER_MOD(moment, zdpl);
double& xxqdp = TINKER_MOD(moment, xxqdp);
double& xyqdp = TINKER_MOD(moment, xyqdp);
double& xzqdp = TINKER_MOD(moment, xzqdp);
double& yxqdp = TINKER_MOD(moment, yxqdp);
double& yyqdp = TINKER_MOD(moment, yyqdp);
double& yzqdp = TINKER_MOD(moment, yzqdp);
double& zxqdp = TINKER_MOD(moment, zxqdp);
double& zyqdp = TINKER_MOD(moment, zyqdp);
double& zzqdp = TINKER_MOD(moment, zzqdp);
#endif

} TINKER_NAMESPACE_END

#endif
