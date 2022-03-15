#pragma once
#include "macro.h"
#include "tool/energy_buffer.h"

namespace tinker {
TINKER_EXTERN int nbond;
TINKER_EXTERN int (*ibnd)[2];
TINKER_EXTERN real* bk;
TINKER_EXTERN real* bl;

TINKER_EXTERN energy_buffer eb;
TINKER_EXTERN virial_buffer vir_eb;
TINKER_EXTERN grad_prec* debx;
TINKER_EXTERN grad_prec* deby;
TINKER_EXTERN grad_prec* debz;
TINKER_EXTERN energy_prec energy_eb;
TINKER_EXTERN virial_prec virial_eb[9];
}
