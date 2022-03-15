#pragma once
#include "elec.h"
#include "mod.chgpen.h"
#include "mod.mplpot.h"

namespace tinker {
void empole_chgpen_data(rc_op op);
void empole_chgpen(int vers);

void empole_chgpen_nonewald(int vers, int use_cf);
void empole_chgpen_ewald(int vers, int use_cf);
void empole_chgpen_ewald_real_self(int vers, int use_cf);
void empole_chgpen_ewald_recip(int vers, int use_cf);

void empole_chgpen_nonewald_acc(int vers, int use_cf);
void empole_chgpen_ewald_recip_acc(int vers, int use_cf);
void empole_chgpen_ewald_real_self_acc(int vers, int use_cf);
void empole_chgpen_nonewald_cu(int vers, int use_cf);
void empole_chgpen_ewald_real_self_cu(int vers, int use_cf);
}
