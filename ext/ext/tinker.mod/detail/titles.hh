#ifndef TINKER_MOD_TITLES_HH_
#define TINKER_MOD_TITLES_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace titles {
extern int& ltitle;
extern char (&title)[240];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(titles, ltitle);
extern "C" char m_tinker_mod(titles, title)[240];

int& ltitle = m_tinker_mod(titles, ltitle);
char (&title)[240] = m_tinker_mod(titles, title);
#endif

} TINKER_NAMESPACE_END

#endif
