#ifndef TINKER_MOD_KIPROP_HH_
#define TINKER_MOD_KIPROP_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kiprop {
const int maxndi = 500;
extern double (&dcon)[maxndi];
extern double (&tdi)[maxndi];
extern char (&kdi)[maxndi][16];

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(kiprop, dcon)[maxndi];
extern "C" double m_tinker_mod(kiprop, tdi)[maxndi];
extern "C" char m_tinker_mod(kiprop, kdi)[maxndi][16];

double (&dcon)[maxndi] = m_tinker_mod(kiprop, dcon);
double (&tdi)[maxndi] = m_tinker_mod(kiprop, tdi);
char (&kdi)[maxndi][16] = m_tinker_mod(kiprop, kdi);
#endif

} TINKER_NAMESPACE_END

#endif
