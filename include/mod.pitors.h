#pragma once
#include "macro.h"
#include "tool/energy_buffer.h"


namespace tinker {
TINKER_EXTERN int npitors;
TINKER_EXTERN int (*ipit)[6];
TINKER_EXTERN real* kpit;


TINKER_EXTERN energy_buffer ept;
TINKER_EXTERN virial_buffer vir_ept;
TINKER_EXTERN grad_prec* deptx;
TINKER_EXTERN grad_prec* depty;
TINKER_EXTERN grad_prec* deptz;
TINKER_EXTERN energy_prec energy_ept;
TINKER_EXTERN virial_prec virial_ept[9];
}
