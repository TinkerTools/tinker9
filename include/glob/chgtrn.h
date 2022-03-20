#pragma once
#include "energybuffer.h"

namespace tinker {
TINKER_EXTERN real* chgct;
TINKER_EXTERN real* dmpct;

TINKER_EXTERN count_buffer nct;
TINKER_EXTERN energy_buffer ect;
TINKER_EXTERN virial_buffer vir_ect;
TINKER_EXTERN grad_prec* dectx;
TINKER_EXTERN grad_prec* decty;
TINKER_EXTERN grad_prec* dectz;
TINKER_EXTERN energy_prec energy_ect;
TINKER_EXTERN virial_prec virial_ect[9];
}
