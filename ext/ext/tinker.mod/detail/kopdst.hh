#ifndef TINKER_MOD_KOPDST_HH_
#define TINKER_MOD_KOPDST_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kopdst {
const int maxnopd = 500;
extern double (&opds)[maxnopd];
extern char (&kopd)[maxnopd][16];

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(kopdst, opds)[maxnopd];
extern "C" char m_tinker_mod(kopdst, kopd)[maxnopd][16];

double (&opds)[maxnopd] = m_tinker_mod(kopdst, opds);
char (&kopd)[maxnopd][16] = m_tinker_mod(kopdst, kopd);
#endif

} TINKER_NAMESPACE_END

#endif
