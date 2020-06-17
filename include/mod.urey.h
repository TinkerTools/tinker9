#pragma once
#include "macro.h"
#include "tool/energy_buffer.h"


namespace tinker {
TINKER_EXTERN int nurey;
TINKER_EXTERN int (*iury)[3];
TINKER_EXTERN real* uk;
TINKER_EXTERN real* ul;


TINKER_EXTERN energy_buffer eub;
TINKER_EXTERN virial_buffer vir_eub;
TINKER_EXTERN grad_prec* deubx;
TINKER_EXTERN grad_prec* deuby;
TINKER_EXTERN grad_prec* deubz;
TINKER_EXTERN energy_prec energy_eub;
TINKER_EXTERN virial_prec virial_eub[9];
}
