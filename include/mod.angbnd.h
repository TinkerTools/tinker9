#pragma once
#include "macro.h"
#include "tool/energybuffer.h"

namespace tinker {
TINKER_EXTERN int nangle;
TINKER_EXTERN int (*iang)[4];
TINKER_EXTERN real* ak;
TINKER_EXTERN real* anat;
TINKER_EXTERN real* afld;

TINKER_EXTERN energy_buffer ea;
TINKER_EXTERN virial_buffer vir_ea;
TINKER_EXTERN grad_prec* deax;
TINKER_EXTERN grad_prec* deay;
TINKER_EXTERN grad_prec* deaz;
TINKER_EXTERN energy_prec energy_ea;
TINKER_EXTERN virial_prec virial_ea[9];
}
