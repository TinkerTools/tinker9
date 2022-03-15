#pragma once
#include "macro.h"

namespace tinker {
enum class ebond_t
{
   harmonic,
   morse
};
TINKER_EXTERN ebond_t bndtyp;
TINKER_EXTERN real bndunit;
TINKER_EXTERN real cbnd;
TINKER_EXTERN real qbnd;
}
