#pragma once
#include "elec.h"
#include "energy_buffer.h"

TINKER_NAMESPACE_BEGIN
TINKER_EXTERN elec_t empole_electyp;

TINKER_EXTERN real m2scale, m3scale, m4scale, m5scale;

TINKER_EXTERN int nmexclude_;
TINKER_EXTERN device_pointer<int, 2> mexclude_;
TINKER_EXTERN device_pointer<real> mexclude_scale_;

TINKER_EXTERN count_buffer nem;
TINKER_EXTERN energy_buffer em;
TINKER_EXTERN virial_buffer vir_em;

void empole_data(rc_op op);


void empole(int vers);
void empole_nonewald(int vers);
void empole_ewald(int vers);
void empole_ewald_real_self(int vers);
void empole_ewald_recip(int vers);

void empole_nonewald_acc(int vers);
void empole_ewald_recip_acc(int vers);
void empole_ewald_real_self_acc(int vers);
void empole_nonewald_cu(int vers);
void empole_ewald_real_self_cu(int vers);
TINKER_NAMESPACE_END
