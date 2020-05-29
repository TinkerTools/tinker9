#pragma once
#include "macro.h"


namespace tinker {
enum class eangle_t : int
{
   in_plane,
   harmonic,
   linear,
   fourier
};
TINKER_EXTERN eangle_t* angtyp;
TINKER_EXTERN real angunit;
TINKER_EXTERN real stbnunit;
TINKER_EXTERN real cang;
TINKER_EXTERN real qang;
TINKER_EXTERN real pang;
TINKER_EXTERN real sang;
}
