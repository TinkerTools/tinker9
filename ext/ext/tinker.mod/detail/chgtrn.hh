#ifndef TINKER_MOD_CHGTRN_HH_
#define TINKER_MOD_CHGTRN_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace chgtrn {
extern int& nct;
extern double*& chgct;
extern double*& dmpct;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(chgtrn, nct);
extern "C" double* m_tinker_mod(chgtrn, chgct);
extern "C" double* m_tinker_mod(chgtrn, dmpct);

int& nct = m_tinker_mod(chgtrn, nct);
double*& chgct = m_tinker_mod(chgtrn, chgct);
double*& dmpct = m_tinker_mod(chgtrn, dmpct);
#endif

} TINKER_NAMESPACE_END

#endif
