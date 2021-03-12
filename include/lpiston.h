#pragma once
#include "mdprec.h"
#include "mod.freeze.h"
#include "tool/energy_buffer.h"


namespace tinker {
extern double lp_alpha;
extern double lp_molpres;
extern double* lp_molpres_buf;
// val = alpha M V**2 -Tr(Xi)
void lp_molpressure(double alpha, double& val);
void lp_molpressure_acc(double alpha, double& val);
void lp_molpressure_cu(double alpha, double& val);


void vv_lpiston_npt(int istep, time_prec dt_ps);
void vv_lpiston_hc_acc(int, time_prec);


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
}
