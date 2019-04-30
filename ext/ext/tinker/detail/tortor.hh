#ifndef TINKER_MOD_TORTOR_HH_
#define TINKER_MOD_TORTOR_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace tortor {
extern int& ntortor;
extern int*& itt;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(tortor, ntortor);
extern "C" int* TINKER_MOD(tortor, itt);

int& ntortor = TINKER_MOD(tortor, ntortor);
int*& itt = TINKER_MOD(tortor, itt);
#endif
} TINKER_NAMESPACE_END

#endif
