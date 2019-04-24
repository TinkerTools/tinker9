#ifndef TINKER_MOD_RXNFLD_HH_
#define TINKER_MOD_RXNFLD_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace rxnfld {
extern int (&ijk)[6][6][6];
extern double (&b1)[13][40];
extern double (&b2)[13][40];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(rxnfld, ijk)[6][6][6];
extern "C" double m_tinker_mod(rxnfld, b1)[13][40];
extern "C" double m_tinker_mod(rxnfld, b2)[13][40];

int (&ijk)[6][6][6] = m_tinker_mod(rxnfld, ijk);
double (&b1)[13][40] = m_tinker_mod(rxnfld, b1);
double (&b2)[13][40] = m_tinker_mod(rxnfld, b2);
#endif

} TINKER_NAMESPACE_END

#endif
