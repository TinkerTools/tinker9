#pragma once
#include "macro.h"
#include "tool/energybuffer.h"

namespace tinker {
TINKER_EXTERN int (*iitors)[4];
TINKER_EXTERN real (*itors1)[4];
TINKER_EXTERN real (*itors2)[4];
TINKER_EXTERN real (*itors3)[4];

TINKER_EXTERN int nitors;
TINKER_EXTERN energy_buffer eit;
TINKER_EXTERN virial_buffer vir_eit;
TINKER_EXTERN grad_prec* deitx;
TINKER_EXTERN grad_prec* deity;
TINKER_EXTERN grad_prec* deitz;
TINKER_EXTERN energy_prec energy_eit;
TINKER_EXTERN virial_prec virial_eit[9];
}
