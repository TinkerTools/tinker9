#ifndef TINKER_MOD_ANGBND_HH_
#define TINKER_MOD_ANGBND_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace angbnd {
extern int& nangle;
extern int*& iang;
extern double*& ak;
extern double*& anat;
extern double*& afld;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(angbnd, nangle);
extern "C" int* m_tinker_mod(angbnd, iang);
extern "C" double* m_tinker_mod(angbnd, ak);
extern "C" double* m_tinker_mod(angbnd, anat);
extern "C" double* m_tinker_mod(angbnd, afld);

int& nangle = m_tinker_mod(angbnd, nangle);
int*& iang = m_tinker_mod(angbnd, iang);
double*& ak = m_tinker_mod(angbnd, ak);
double*& anat = m_tinker_mod(angbnd, anat);
double*& afld = m_tinker_mod(angbnd, afld);
#endif

} TINKER_NAMESPACE_END

#endif
