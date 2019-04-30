#ifndef TINKER_MOD_TITLES_HH_
#define TINKER_MOD_TITLES_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace titles {
extern int& ltitle;
extern char (&title)[240];

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(titles, ltitle);
extern "C" char TINKER_MOD(titles, title)[240];

int& ltitle = TINKER_MOD(titles, ltitle);
char (&title)[240] = TINKER_MOD(titles, title);
#endif
} TINKER_NAMESPACE_END

#endif
