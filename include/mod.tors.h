#pragma once
#include "macro.h"
#include "tool/energybuffer.h"

namespace tinker {
TINKER_EXTERN int ntors;
TINKER_EXTERN int (*itors)[4];
TINKER_EXTERN real (*tors1)[4];
TINKER_EXTERN real (*tors2)[4];
TINKER_EXTERN real (*tors3)[4];
TINKER_EXTERN real (*tors4)[4];
TINKER_EXTERN real (*tors5)[4];
TINKER_EXTERN real (*tors6)[4];

TINKER_EXTERN energy_buffer et;
TINKER_EXTERN virial_buffer vir_et;
TINKER_EXTERN grad_prec* detx;
TINKER_EXTERN grad_prec* dety;
TINKER_EXTERN grad_prec* detz;
TINKER_EXTERN energy_prec energy_et;
TINKER_EXTERN virial_prec virial_et[9];
}
