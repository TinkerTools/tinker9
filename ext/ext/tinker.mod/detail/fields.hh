#ifndef TINKER_MOD_FIELDS_HH_
#define TINKER_MOD_FIELDS_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace fields {
const int maxbio = 10000;
extern int (&biotyp)[maxbio];
extern char (&forcefield)[20];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(fields, biotyp)[maxbio];
extern "C" char m_tinker_mod(fields, forcefield)[20];

int (&biotyp)[maxbio] = m_tinker_mod(fields, biotyp);
char (&forcefield)[20] = m_tinker_mod(fields, forcefield);
#endif

} TINKER_NAMESPACE_END

#endif
