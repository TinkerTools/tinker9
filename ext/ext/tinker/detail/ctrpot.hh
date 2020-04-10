#pragma once

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace ctrpot {
extern char (&ctrntyp)[8];

#ifdef TINKER_MOD_CPP_
extern "C" char TINKER_MOD(ctrpot, ctrntyp)[8];

char (&ctrntyp)[8] = TINKER_MOD(ctrpot, ctrntyp);
#endif
} TINKER_NAMESPACE_END
