#pragma once
#include "macro.h"
#include "tool/energybuffer.h"

namespace tinker {
TINKER_EXTERN int ntortor;
TINKER_EXTERN int (*itt)[3];

//====================================================================//

TINKER_EXTERN int* chkttor_ia_;

TINKER_EXTERN energy_buffer ett;
TINKER_EXTERN virial_buffer vir_ett;
TINKER_EXTERN grad_prec* dettx;
TINKER_EXTERN grad_prec* detty;
TINKER_EXTERN grad_prec* dettz;
TINKER_EXTERN energy_prec energy_ett;
TINKER_EXTERN virial_prec virial_ett[9];
}
