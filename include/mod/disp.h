#pragma once
#include "ff/energybuffer.h"

namespace tinker {
TINKER_EXTERN real csixpr;
TINKER_EXTERN real* csix;
TINKER_EXTERN real* adisp;

TINKER_EXTERN int ndspexclude;
TINKER_EXTERN int (*dspexclude)[2];
TINKER_EXTERN real* dspexclude_scale;

TINKER_EXTERN CountBuffer ndisp;
TINKER_EXTERN EnergyBuffer edsp;
TINKER_EXTERN VirialBuffer vir_edsp;
TINKER_EXTERN grad_prec* dedspx;
TINKER_EXTERN grad_prec* dedspy;
TINKER_EXTERN grad_prec* dedspz;
TINKER_EXTERN energy_prec energy_edsp;
TINKER_EXTERN virial_prec virial_edsp[9];

TINKER_EXTERN real dsp2scale;
TINKER_EXTERN real dsp3scale;
TINKER_EXTERN real dsp4scale;
TINKER_EXTERN real dsp5scale;

TINKER_EXTERN real elrc_vol_dsp;
TINKER_EXTERN real vlrc_vol_dsp;
}
