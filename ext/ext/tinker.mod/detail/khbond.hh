#ifndef TINKER_MOD_KHBOND_HH_
#define TINKER_MOD_KHBOND_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace khbond {
const int maxnhb = 500;
extern double (&radhb)[maxnhb];
extern double (&epshb)[maxnhb];
extern char (&khb)[maxnhb][8];

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(khbond, radhb)[maxnhb];
extern "C" double m_tinker_mod(khbond, epshb)[maxnhb];
extern "C" char m_tinker_mod(khbond, khb)[maxnhb][8];

double (&radhb)[maxnhb] = m_tinker_mod(khbond, radhb);
double (&epshb)[maxnhb] = m_tinker_mod(khbond, epshb);
char (&khb)[maxnhb][8] = m_tinker_mod(khbond, khb);
#endif

} TINKER_NAMESPACE_END

#endif
