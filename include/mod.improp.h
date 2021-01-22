#pragma once
#include "macro.h"
#include "tool/energy_buffer.h"


namespace tinker {
TINKER_EXTERN int (*iiprop)[4];
TINKER_EXTERN real* kprop;
TINKER_EXTERN real* vprop;


TINKER_EXTERN int niprop;
TINKER_EXTERN energy_buffer eid;
TINKER_EXTERN virial_buffer vir_eid;
TINKER_EXTERN grad_prec* deidx;
TINKER_EXTERN grad_prec* deidy;
TINKER_EXTERN grad_prec* deidz;
TINKER_EXTERN energy_prec energy_eid;
TINKER_EXTERN virial_prec virial_eid[9];
}
