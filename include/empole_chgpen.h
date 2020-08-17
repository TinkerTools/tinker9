#pragma once
#include "elec.h"
#include "mod.mplpot.h"
#include "mod.chgpen.h"



namespace tinker {
void empole_chgpen_data(rc_op op);
void empole_chgpen(int vers);


void empole_chgpen_nonewald(int vers);
void empole_chgpen_ewald(int vers);
void empole_chgpen_ewald_real_self(int vers);
void empole_chgpen_ewald_recip(int vers);


void empole_chgpen_nonewald_acc(int vers);
void empole_chgpen_ewald_recip_acc(int vers);
void empole_chgpen_ewald_real_self_acc(int vers);
void empole_chgpen_nonewald_cu(int vers);
void empole_chgpen_ewald_real_self_cu(int vers);
}
