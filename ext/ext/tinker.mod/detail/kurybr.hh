#ifndef TINKER_MOD_KURYBR_HH_
#define TINKER_MOD_KURYBR_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kurybr {
const int maxnu = 2000;
extern double (&ucon)[maxnu];
extern double (&dst13)[maxnu];
extern char (&ku)[maxnu][12];

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(kurybr, ucon)[maxnu];
extern "C" double m_tinker_mod(kurybr, dst13)[maxnu];
extern "C" char m_tinker_mod(kurybr, ku)[maxnu][12];

double (&ucon)[maxnu] = m_tinker_mod(kurybr, ucon);
double (&dst13)[maxnu] = m_tinker_mod(kurybr, dst13);
char (&ku)[maxnu][12] = m_tinker_mod(kurybr, ku);
#endif

} TINKER_NAMESPACE_END

#endif
