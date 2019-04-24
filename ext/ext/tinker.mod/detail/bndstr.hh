#ifndef TINKER_MOD_BNDSTR_HH_
#define TINKER_MOD_BNDSTR_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace bndstr {
extern int& nbond;
extern int*& ibnd;
extern double*& bk;
extern double*& bl;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(bndstr, nbond);
extern "C" int* m_tinker_mod(bndstr, ibnd);
extern "C" double* m_tinker_mod(bndstr, bk);
extern "C" double* m_tinker_mod(bndstr, bl);

int& nbond = m_tinker_mod(bndstr, nbond);
int*& ibnd = m_tinker_mod(bndstr, ibnd);
double*& bk = m_tinker_mod(bndstr, bk);
double*& bl = m_tinker_mod(bndstr, bl);
#endif

} TINKER_NAMESPACE_END

#endif
