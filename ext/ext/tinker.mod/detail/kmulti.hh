#ifndef TINKER_MOD_KMULTI_HH_
#define TINKER_MOD_KMULTI_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kmulti {
const int maxnmp = 2000;
extern double (&multip)[maxnmp][13];
extern char (&mpaxis)[maxnmp][8];
extern char (&kmp)[maxnmp][16];

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(kmulti, multip)[maxnmp][13];
extern "C" char m_tinker_mod(kmulti, mpaxis)[maxnmp][8];
extern "C" char m_tinker_mod(kmulti, kmp)[maxnmp][16];

double (&multip)[maxnmp][13] = m_tinker_mod(kmulti, multip);
char (&mpaxis)[maxnmp][8] = m_tinker_mod(kmulti, mpaxis);
char (&kmp)[maxnmp][16] = m_tinker_mod(kmulti, kmp);
#endif

} TINKER_NAMESPACE_END

#endif
