#ifndef TINKER_MOD_KSTTOR_HH_
#define TINKER_MOD_KSTTOR_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace ksttor {
const int maxnbt = 500;
extern double (&btcon)[maxnbt][9];
extern char (&kbt)[maxnbt][16];

#ifdef TINKER_MOD_CPP_
extern "C" double TINKER_MOD(ksttor, btcon)[maxnbt][9];
extern "C" char TINKER_MOD(ksttor, kbt)[maxnbt][16];

double (&btcon)[maxnbt][9] = TINKER_MOD(ksttor, btcon);
char (&kbt)[maxnbt][16] = TINKER_MOD(ksttor, kbt);
#endif

} TINKER_NAMESPACE_END

#endif
