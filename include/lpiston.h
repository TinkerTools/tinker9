#pragma once
#include "md.h"
#include "mod.freeze.h"
#include "tool/energybuffer.h"

namespace tinker {
extern double lp_rats1;
extern double lp_eksum;
extern double lp_ekin[3][3];
extern double lp_vir[9];
extern virial_buffer lp_vir_buf;

extern double vbar_matrix[3][3];

// (ax,ay,az) = mat (ax,ay,az) or (ax,ay,az) = transpose(mat) (ax,ay,az)
void lp_matvec(int len, char transpose, double mat[3][3], pos_prec* ax, pos_prec* ay, pos_prec* az);

void lp_atom_kinetic();
void lp_mol_kinetic();
void lp_virial(bool molP);

void propagate_pos_raxbv(double a, double b);
void propagate_pos_raxbv_aniso(double a[3][3], double b[3][3]);
void propagate_pos_axbv_aniso(double a[3][3], double b[3][3]);
void lp_propagate_mol_vel(vel_prec scal);
void lp_propagate_mol_vel_aniso(vel_prec scal[3][3]);
void lp_center_of_mass(const pos_prec* atomx, const pos_prec* atomy, const pos_prec* atomz,
   pos_prec* molx, pos_prec* moly, pos_prec* molz);

void lprat(time_prec dt, const pos_prec* xold, const pos_prec* yold, const pos_prec* zold);

void vv_lpiston_init();
void vv_lpiston_destory();
void vv_lpiston_npt(int istep, time_prec dt_ps);
}
