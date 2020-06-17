#pragma once
#include "macro.h"
#include "tool/energy_buffer.h"


namespace tinker {
// valence term gradients also work as total gradients
TINKER_EXTERN grad_prec* gx;
TINKER_EXTERN grad_prec* gy;
TINKER_EXTERN grad_prec* gz;


// AMOEBA vdw, HIPPO dispersion, HIPPO repulsion
TINKER_EXTERN grad_prec* gx_vdw;
TINKER_EXTERN grad_prec* gy_vdw;
TINKER_EXTERN grad_prec* gz_vdw;


// AMOEBA charge, AMOEBA mpole, AMOEBA polar
// HIPPO mpole, HIPPO polar, HIPPO charge transfer
TINKER_EXTERN grad_prec* gx_elec;
TINKER_EXTERN grad_prec* gy_elec;
TINKER_EXTERN grad_prec* gz_elec;


TINKER_EXTERN energy_buffer eng_buf;
TINKER_EXTERN energy_buffer eng_buf_vdw;
TINKER_EXTERN energy_buffer eng_buf_elec;


TINKER_EXTERN virial_buffer vir_buf;
TINKER_EXTERN virial_buffer vir_buf_vdw;
TINKER_EXTERN virial_buffer vir_buf_elec;


TINKER_EXTERN energy_prec esum; // total potential energy
TINKER_EXTERN energy_prec energy_valence;
TINKER_EXTERN energy_prec energy_vdw;
TINKER_EXTERN energy_prec energy_elec;


TINKER_EXTERN virial_prec vir[9]; // total virial tensor
TINKER_EXTERN virial_prec virial_valence[9];
TINKER_EXTERN virial_prec virial_vdw[9];
TINKER_EXTERN virial_prec virial_elec[9];
}
