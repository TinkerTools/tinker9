#pragma once
#include "macro.h"

namespace tinker {
TINKER_EXTERN real m2scale;
TINKER_EXTERN real m3scale;
TINKER_EXTERN real m4scale;
TINKER_EXTERN real m5scale;

//====================================================================//

TINKER_EXTERN int nmexclude;
TINKER_EXTERN int (*mexclude)[2];
TINKER_EXTERN real* mexclude_scale;
}
