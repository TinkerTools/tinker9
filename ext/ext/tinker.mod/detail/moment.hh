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
extern "C" double m_tinker_mod(moment, netchg);
extern "C" double m_tinker_mod(moment, netdpl);
extern "C" double m_tinker_mod(moment, netqdp)[3];
extern "C" double m_tinker_mod(moment, xdpl);
extern "C" double m_tinker_mod(moment, ydpl);
extern "C" double m_tinker_mod(moment, zdpl);
extern "C" double m_tinker_mod(moment, xxqdp);
extern "C" double m_tinker_mod(moment, xyqdp);
extern "C" double m_tinker_mod(moment, xzqdp);
extern "C" double m_tinker_mod(moment, yxqdp);
extern "C" double m_tinker_mod(moment, yyqdp);
extern "C" double m_tinker_mod(moment, yzqdp);
extern "C" double m_tinker_mod(moment, zxqdp);
extern "C" double m_tinker_mod(moment, zyqdp);
extern "C" double m_tinker_mod(moment, zzqdp);

double& netchg = m_tinker_mod(moment, netchg);
double& netdpl = m_tinker_mod(moment, netdpl);
double (&netqdp)[3] = m_tinker_mod(moment, netqdp);
double& xdpl = m_tinker_mod(moment, xdpl);
double& ydpl = m_tinker_mod(moment, ydpl);
double& zdpl = m_tinker_mod(moment, zdpl);
double& xxqdp = m_tinker_mod(moment, xxqdp);
double& xyqdp = m_tinker_mod(moment, xyqdp);
double& xzqdp = m_tinker_mod(moment, xzqdp);
double& yxqdp = m_tinker_mod(moment, yxqdp);
double& yyqdp = m_tinker_mod(moment, yyqdp);
double& yzqdp = m_tinker_mod(moment, yzqdp);
double& zxqdp = m_tinker_mod(moment, zxqdp);
double& zyqdp = m_tinker_mod(moment, zyqdp);
double& zzqdp = m_tinker_mod(moment, zzqdp);
#endif

} TINKER_NAMESPACE_END

#endif
