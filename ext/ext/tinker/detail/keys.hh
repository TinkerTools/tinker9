#ifndef TINKER_MOD_KEYS_HH_
#define TINKER_MOD_KEYS_HH_

#include "util_macro.h"

TINKER_NAMESPACE_BEGIN namespace keys {
const int maxkey = 25000;
extern int& nkey;
extern char (&keyline)[maxkey][240];

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(keys, nkey);
extern "C" char TINKER_MOD(keys, keyline)[maxkey][240];

int& nkey = TINKER_MOD(keys, nkey);
char (&keyline)[maxkey][240] = TINKER_MOD(keys, keyline);
#endif
} TINKER_NAMESPACE_END

#endif
