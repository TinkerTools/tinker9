#ifndef TINKER_MOD_OPDIST_HH_
#define TINKER_MOD_OPDIST_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace opdist {
extern int& nopdist;
extern int*& iopd;
extern double*& opdk;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(opdist, nopdist);
extern "C" int* m_tinker_mod(opdist, iopd);
extern "C" double* m_tinker_mod(opdist, opdk);

int& nopdist = m_tinker_mod(opdist, nopdist);
int*& iopd = m_tinker_mod(opdist, iopd);
double*& opdk = m_tinker_mod(opdist, opdk);
#endif

} TINKER_NAMESPACE_END

#endif
