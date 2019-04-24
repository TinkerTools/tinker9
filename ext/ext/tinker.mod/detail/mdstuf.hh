#ifndef TINKER_MOD_MDSTUF_HH_
#define TINKER_MOD_MDSTUF_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace mdstuf {
extern int& nfree;
extern int& irest;
extern int& bmnmix;
extern double& arespa;
extern int& dorest;
extern char (&integrate)[11];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(mdstuf, nfree);
extern "C" int m_tinker_mod(mdstuf, irest);
extern "C" int m_tinker_mod(mdstuf, bmnmix);
extern "C" double m_tinker_mod(mdstuf, arespa);
extern "C" int m_tinker_mod(mdstuf, dorest);
extern "C" char m_tinker_mod(mdstuf, integrate)[11];

int& nfree = m_tinker_mod(mdstuf, nfree);
int& irest = m_tinker_mod(mdstuf, irest);
int& bmnmix = m_tinker_mod(mdstuf, bmnmix);
double& arespa = m_tinker_mod(mdstuf, arespa);
int& dorest = m_tinker_mod(mdstuf, dorest);
char (&integrate)[11] = m_tinker_mod(mdstuf, integrate);
#endif

} TINKER_NAMESPACE_END

#endif
