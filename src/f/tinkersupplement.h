#pragma once
#include "tool/macro.h"
#include <tinker/routines.h>

#define TINKER_SUPPL_DECL

#ifdef __cplusplus
TINKER_DECL_EXTN("C")
{
#endif

#include "suppl.f"

#ifdef __cplusplus
}
#endif
