#pragma once
#include "macro.h"
#include "tool/energybuffer.h"

namespace tinker {
/**
 * \ingroup charge
 * \brief Magnitude of the partial charges (e-) of each **atom**.
 * \note Unlike Tinker, where only non-zero charges will be stored, this
 * array also includes zero charges.
 */
TINKER_EXTERN real* pchg;

TINKER_EXTERN count_buffer nec;
TINKER_EXTERN energy_buffer ec;
TINKER_EXTERN virial_buffer vir_ec;
TINKER_EXTERN grad_prec* decx;
TINKER_EXTERN grad_prec* decy;
TINKER_EXTERN grad_prec* decz;
TINKER_EXTERN energy_prec energy_ec;
TINKER_EXTERN virial_prec virial_ec[9];
}
