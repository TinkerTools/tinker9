#pragma once
#include "ff/energybuffer.h"

// charge chgpot
namespace tinker {
/// \ingroup charge
/// \brief Magnitude of the partial charges (e-) of each **atom**.
/// \note Unlike Tinker, where only non-zero charges will be stored, this
/// array also includes zero charges.
TINKER_EXTERN real* pchg;

TINKER_EXTERN count_buffer nec;
TINKER_EXTERN energy_buffer ec;
TINKER_EXTERN virial_buffer vir_ec;
TINKER_EXTERN grad_prec* decx;
TINKER_EXTERN grad_prec* decy;
TINKER_EXTERN grad_prec* decz;
TINKER_EXTERN energy_prec energy_ec;
TINKER_EXTERN virial_prec virial_ec[9];

TINKER_EXTERN real ebuffer;
TINKER_EXTERN real c2scale, c3scale, c4scale, c5scale;
TINKER_EXTERN int ncexclude;
TINKER_EXTERN int (*cexclude)[2];
TINKER_EXTERN real* cexclude_scale;
}

// chglj
namespace tinker {
TINKER_EXTERN int ncvexclude;
TINKER_EXTERN int (*cvexclude)[2];
TINKER_EXTERN real (*cvexclude_scale)[2];

TINKER_EXTERN bool vdwpr_in_use;
}
