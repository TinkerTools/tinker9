#pragma once
#include "mdcalc.h"
#include "mdegv.h"
#include "mdprec.h"
#include "tool/energy_buffer.h"


namespace tinker {
// clang-format off
#define TINKER_COUNT_BUFFERS \
   nev, \
   nec, nem, nep, \
   nct
#define TINKER_ENERGY_BUFFERS \
   eb, ea, eba, eub, \
   eopb, et, ept, ett, \
   eg, \
   \
   ev, \
   ec, em, ep, \
   ect, \
   \
   eng_buf
#define TINKER_VIRIAL_BUFFERS \
   vir_eb, vir_ea, vir_eba, vir_eub, \
   vir_eopb, vir_et, vir_ept, vir_ett, \
   vir_eg, \
   \
   vir_ev, \
   vir_ec, vir_em, vir_ep, vir_trq, \
   vir_ect, \
   \
   vir_buf
#define TINKER_GRADIENTS \
   debx, deby, debz, \
   deax, deay, deaz, \
   debax, debay, debaz, \
   deubx, deuby, deubz, \
   deopbx, deopby, deopbz, \
   detx, dety, detz, \
   deptx, depty, deptz, \
   dettx, detty, dettz, \
   degx, degy, degz, \
   \
   devx, devy, devz, \
   decx, decy, decz, \
   demx, demy, demz, \
   depx, depy, depz, \
   \
   dectx, decty, dectz, \
   \
   gx, gy, gz


#define TINKER_ENERGY_VARIABLES  \
   energy_eb, energy_ea, energy_eba, energy_eub, \
   energy_eopb, energy_et, energy_ept, energy_ett, \
   energy_eg, \
   \
   energy_ev, \
   energy_ec, energy_em, energy_ep, \
   energy_ect, \
   \
   esum, \
   energy_valence, energy_elec
#define TINKER_VIRIAL_TENSORS \
   virial_eb, virial_ea, virial_eba, virial_eub, \
   virial_eopb, virial_et, virial_ept, virial_ett, \
   virial_eg, \
   \
   virial_ev, \
   \
   virial_ec, virial_em, virial_ep,\
   \
   vir, \
   virial_valence, virial_elec
// clang-format on


TINKER_EXTERN count_buffer TINKER_COUNT_BUFFERS;
TINKER_EXTERN energy_buffer TINKER_ENERGY_BUFFERS;
TINKER_EXTERN virial_buffer TINKER_VIRIAL_BUFFERS;
using _detail_GradPrecPointer = grad_prec*;
TINKER_EXTERN _detail_GradPrecPointer TINKER_GRADIENTS;
TINKER_EXTERN energy_prec TINKER_ENERGY_VARIABLES;
using _detail_VirialTensor = virial_prec[9];
TINKER_EXTERN _detail_VirialTensor TINKER_VIRIAL_TENSORS;
}
