#ifndef TINKER_MOD_KANANG_HH_
#define TINKER_MOD_KANANG_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kanang {
extern double*& anan;

#ifdef TINKER_MOD_CPP_
extern "C" double* m_tinker_mod(kanang, anan);

double*& anan = m_tinker_mod(kanang, anan);
#endif

} TINKER_NAMESPACE_END

#endif
