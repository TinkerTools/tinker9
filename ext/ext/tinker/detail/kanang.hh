#ifndef TINKER_MOD_KANANG_HH_
#define TINKER_MOD_KANANG_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace kanang {
extern double*& anan;

#ifdef TINKER_MOD_CPP_
extern "C" double* TINKER_MOD(kanang, anan);

double*& anan = TINKER_MOD(kanang, anan);
#endif
} TINKER_NAMESPACE_END

#endif
