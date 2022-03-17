#pragma once
#include "macro.h"
#include "tool/energybuffer.h"

namespace tinker {
TINKER_EXTERN int nstrtor;
TINKER_EXTERN int (*ist)[4];
TINKER_EXTERN real (*kst)[9];

TINKER_EXTERN energy_buffer ebt;
TINKER_EXTERN virial_buffer vir_ebt;
TINKER_EXTERN grad_prec* debtx;
TINKER_EXTERN grad_prec* debty;
TINKER_EXTERN grad_prec* debtz;
TINKER_EXTERN energy_prec energy_ebt;
TINKER_EXTERN virial_prec virial_ebt[9];
}
