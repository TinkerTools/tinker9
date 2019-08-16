#ifndef TINKER_MOD_RXNFLD_HH_
#define TINKER_MOD_RXNFLD_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace rxnfld {
extern int (&ijk)[6][6][6];
extern double (&b1)[13][40];
extern double (&b2)[13][40];

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(rxnfld, ijk)[6][6][6];
extern "C" double TINKER_MOD(rxnfld, b1)[13][40];
extern "C" double TINKER_MOD(rxnfld, b2)[13][40];

int (&ijk)[6][6][6] = TINKER_MOD(rxnfld, ijk);
double (&b1)[13][40] = TINKER_MOD(rxnfld, b1);
double (&b2)[13][40] = TINKER_MOD(rxnfld, b2);
#endif
} TINKER_NAMESPACE_END

#endif
