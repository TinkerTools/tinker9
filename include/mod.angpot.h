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


enum class eopbend_t
{
   w_d_c,
   allinger
};


TINKER_EXTERN real angunit;
TINKER_EXTERN real stbnunit;
TINKER_EXTERN real opbunit;
TINKER_EXTERN real cang;
TINKER_EXTERN real qang;
TINKER_EXTERN real pang;
TINKER_EXTERN real sang;
TINKER_EXTERN real copb;
TINKER_EXTERN real qopb;
TINKER_EXTERN real popb;
TINKER_EXTERN real sopb;
TINKER_EXTERN eopbend_t opbtyp;
TINKER_EXTERN eangle_t* angtyp;
}
