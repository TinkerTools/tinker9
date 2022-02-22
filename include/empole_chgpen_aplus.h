#pragma once
#include "elec.h"
#include "mod.chgpen.h"
#include "mod.mplpot.h"


namespace tinker {
void empole_chgpen_aplus_data(rc_op op);
void empole_chgpen_aplus(int vers);


void empole_chgpen_aplus_nonewald(int vers, int use_cf);
void empole_chgpen_aplus_ewald(int vers, int use_cf);
void empole_chgpen_aplus_ewald_real_self(int vers, int use_cf);
void empole_chgpen_aplus_ewald_recip(int vers, int use_cf);


void empole_chgpen_aplus_nonewald_acc(int vers, int use_cf);
void empole_chgpen_aplus_ewald_recip_acc(int vers, int use_cf);
void empole_chgpen_aplus_ewald_real_self_acc(int vers, int use_cf);
void empole_chgpen_aplus_nonewald_cu(int vers, int use_cf);
void empole_chgpen_aplus_ewald_real_self_cu(int vers, int use_cf);
}
