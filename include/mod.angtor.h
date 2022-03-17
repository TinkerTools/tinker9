#pragma once
#include "macro.h"
#include "tool/energybuffer.h"

namespace tinker {
TINKER_EXTERN int nangtor;
TINKER_EXTERN int (*iat)[3];
TINKER_EXTERN real (*kant)[6];

TINKER_EXTERN energy_buffer eat;
TINKER_EXTERN virial_buffer vir_eat;
TINKER_EXTERN grad_prec* deatx;
TINKER_EXTERN grad_prec* deaty;
TINKER_EXTERN grad_prec* deatz;
TINKER_EXTERN energy_prec energy_eat;
TINKER_EXTERN virial_prec virial_eat[9];
}
