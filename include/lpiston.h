#pragma once
#include "mdprec.h"
#include "mod.freeze.h"


namespace tinker {
void vv_lpiston_npt(int istep, time_prec dt_ps);


void lprat(time_prec dt, const pos_prec* xold, const pos_prec* yold,
           const pos_prec* zold);
void lprat_acc(time_prec, const pos_prec*, const pos_prec*, const pos_prec*);
void lprat_settle_acc(time_prec, const pos_prec*, const pos_prec*,
                      const pos_prec*);
void lprat_ch_acc(time_prec, const pos_prec*, const pos_prec*, const pos_prec*);
void lprat_methyl_cu(time_prec, const pos_prec*, const pos_prec*,
                     const pos_prec*);


// always compute virial (atomic and molecular virial)
void lprat2(time_prec dt);
void lprat2_acc(time_prec);
void lprat2_settle_acc(time_prec);
void lprat2_ch_acc(time_prec);
void lprat2_methyl_cu(time_prec);


// val = coef * mol_kinetic + mol_vir
void ratcom_kevir(double coef, double atomic_vir, double& val);
void ratcom_kevir_acc(double coef, double atomic_vir, double& val);
void ratcom_kevir_cu(double coef, double atomic_vir, double& val);
}
