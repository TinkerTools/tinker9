#pragma once
#include "macro.h"
#include "tool/energybuffer.h"

namespace tinker {
TINKER_EXTERN int nstrbnd;
TINKER_EXTERN int (*isb)[3];
TINKER_EXTERN real (*sbk)[2];

TINKER_EXTERN energy_buffer eba;
TINKER_EXTERN virial_buffer vir_eba;
TINKER_EXTERN grad_prec* debax;
TINKER_EXTERN grad_prec* debay;
TINKER_EXTERN grad_prec* debaz;
TINKER_EXTERN energy_prec energy_eba;
TINKER_EXTERN virial_prec virial_eba[9];
}
