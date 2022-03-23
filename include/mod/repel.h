#pragma once
#include "ff/energybuffer.h"

namespace tinker {
TINKER_EXTERN real* sizpr;
TINKER_EXTERN real* dmppr;
TINKER_EXTERN real* elepr;

TINKER_EXTERN int nrepexclude;
TINKER_EXTERN int (*repexclude)[2];
TINKER_EXTERN real* repexclude_scale;

TINKER_EXTERN count_buffer nrep;
TINKER_EXTERN energy_buffer er;
TINKER_EXTERN virial_buffer vir_er;
TINKER_EXTERN grad_prec* derx;
TINKER_EXTERN grad_prec* dery;
TINKER_EXTERN grad_prec* derz;
TINKER_EXTERN energy_prec energy_er;
TINKER_EXTERN virial_prec virial_er[9];

TINKER_EXTERN real r2scale;
TINKER_EXTERN real r3scale;
TINKER_EXTERN real r4scale;
TINKER_EXTERN real r5scale;
}
