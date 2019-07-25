#ifndef TINKER_MOD_CHGTRN_HH_
#define TINKER_MOD_CHGTRN_HH_

#include "util_macro.h"

TINKER_NAMESPACE_BEGIN namespace chgtrn {
extern int& nct;
extern double*& chgct;
extern double*& dmpct;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(chgtrn, nct);
extern "C" double* TINKER_MOD(chgtrn, chgct);
extern "C" double* TINKER_MOD(chgtrn, dmpct);

int& nct = TINKER_MOD(chgtrn, nct);
double*& chgct = TINKER_MOD(chgtrn, chgct);
double*& dmpct = TINKER_MOD(chgtrn, dmpct);
#endif
} TINKER_NAMESPACE_END

#endif
