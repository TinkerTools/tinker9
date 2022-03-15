#pragma once
#include "macro.h"
#include "tool/energy_buffer.h"

namespace tinker {
TINKER_EXTERN int nopbend;
TINKER_EXTERN int* iopb;
TINKER_EXTERN real* opbk;

TINKER_EXTERN energy_buffer eopb;
TINKER_EXTERN virial_buffer vir_eopb;
TINKER_EXTERN grad_prec* deopbx;
TINKER_EXTERN grad_prec* deopby;
TINKER_EXTERN grad_prec* deopbz;
TINKER_EXTERN energy_prec energy_eopb;
TINKER_EXTERN virial_prec virial_eopb[9];
}
