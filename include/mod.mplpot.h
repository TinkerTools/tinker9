#pragma once
#include "macro.h"


namespace tinker {
TINKER_EXTERN real m2scale;
TINKER_EXTERN real m3scale;
TINKER_EXTERN real m4scale;
TINKER_EXTERN real m5scale;

enum class chgpen_t
{
   GORDON1,
   GORDON2
};
TINKER_EXTERN chgpen_t pentyp;


//====================================================================//


TINKER_EXTERN int nmexclude;
TINKER_EXTERN int (*mexclude)[2];
TINKER_EXTERN real* mexclude_scale;
}
