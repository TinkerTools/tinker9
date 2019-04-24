#ifndef TINKER_MOD_KEYS_HH_
#define TINKER_MOD_KEYS_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace keys {
const int maxkey = 25000;
extern int& nkey;
extern char (&keyline)[maxkey][240];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(keys, nkey);
extern "C" char m_tinker_mod(keys, keyline)[maxkey][240];

int& nkey = m_tinker_mod(keys, nkey);
char (&keyline)[maxkey][240] = m_tinker_mod(keys, keyline);
#endif

} TINKER_NAMESPACE_END

#endif
