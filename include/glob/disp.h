#pragma once
#include "energybuffer.h"

namespace tinker {
TINKER_EXTERN real csixpr;
TINKER_EXTERN real* csix;
TINKER_EXTERN real* adisp;

TINKER_EXTERN int ndspexclude;
TINKER_EXTERN int (*dspexclude)[2];
TINKER_EXTERN real* dspexclude_scale;

TINKER_EXTERN count_buffer ndisp;
TINKER_EXTERN energy_buffer edsp;
TINKER_EXTERN virial_buffer vir_edsp;
TINKER_EXTERN grad_prec* dedspx;
TINKER_EXTERN grad_prec* dedspy;
TINKER_EXTERN grad_prec* dedspz;
TINKER_EXTERN energy_prec energy_edsp;
TINKER_EXTERN virial_prec virial_edsp[9];
}
