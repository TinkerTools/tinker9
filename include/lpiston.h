#pragma once
#include "mdprec.h"
#include "mod.freeze.h"
#include "tool/energy_buffer.h"


namespace tinker {
extern double lp_alpha;
// extern double lp_mol_eksum;
// extern double lp_mol_trvir;
// extern double* lp_mol_ek_buf;
// extern double* lp_mol_vir_buf;
// extern double* lp_mol_ekvir_buf;


extern double lp_rats1, lp_rats2;


extern double lp_eksum;
extern double lp_ekin[3][3];
extern double lp_vir[9];
extern virial_buffer lp_vir_buf;
void lp_atom_kinetic();
void lp_mol_kinetic();
void lp_atom_virial();
void lp_mol_virial();


void vv_lpiston_init();
void vv_lpiston_destory();


void vv_lpiston_npt(int istep, time_prec dt_ps);
void vv_lpiston_npt_acc(int, time_prec);

void lp_center_of_mass(const pos_prec* atomx, const pos_prec* atomy,
                       const pos_prec* atomz, pos_prec* molx, pos_prec* moly,
                       pos_prec* molz);
void lp_center_of_mass_acc(const pos_prec*, const pos_prec*, const pos_prec*,
                           pos_prec*, pos_prec*, pos_prec*);


void lprat(time_prec dt, const pos_prec* xold, const pos_prec* yold,
           const pos_prec* zold);
void lprat_acc(time_prec, const pos_prec*, const pos_prec*, const pos_prec*);
void lprat_settle_acc(time_prec, const pos_prec*, const pos_prec*,
                      const pos_prec*);
void lprat_ch_acc(time_prec, const pos_prec*, const pos_prec*, const pos_prec*);
void lprat_methyl_cu(time_prec, const pos_prec*, const pos_prec*,
                     const pos_prec*);


void lprat2(time_prec dt);
void lprat2_acc(time_prec);
void lprat2_settle_acc(time_prec);
void lprat2_ch_acc(time_prec);
void lprat2_methyl_cu(time_prec);
}
