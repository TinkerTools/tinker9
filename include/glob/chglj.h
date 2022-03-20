#pragma once
#include "macro.h"

namespace tinker {
TINKER_EXTERN int ncvexclude;
TINKER_EXTERN int (*cvexclude)[2];
TINKER_EXTERN real (*cvexclude_scale)[2];

TINKER_EXTERN bool vdwpr_in_use;
}
