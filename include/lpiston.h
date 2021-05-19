#pragma once
#include "mdprec.h"
#include "mod.freeze.h"
#include "tool/energy_buffer.h"


namespace tinker {
extern double lp_alpha;
extern double lp_rats1;
extern double lp_rats2;
extern double lp_eksum;
extern double lp_ekin[3][3];
extern double lp_vir[9];
extern virial_buffer lp_vir_buf;


void lp_atom_kinetic();
void lp_mol_kinetic();
void lp_atom_virial();
void lp_mol_virial();


double sinh_id(double x);
void lp_center_of_mass(const pos_prec* atomx, const pos_prec* atomy,
                       const pos_prec* atomz, pos_prec* molx, pos_prec* moly,
                       pos_prec* molz);


void lprat(time_prec dt, const pos_prec* xold, const pos_prec* yold,
           const pos_prec* zold);
void lprat2(time_prec dt);


void vv_lpiston_init();
void vv_lpiston_destory();
void vv_lpiston_npt(int istep, time_prec dt_ps);
}
