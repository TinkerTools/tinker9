#pragma once
#include "mod.freeze.h"
#include "tool/rc_man.h"


namespace tinker {
bool use_rattle();
void rattle_data(rc_op);


// val = coef * mol_kinetic + mol_vir
void ratcom_kevir(double coef, double atomic_vir, double& val);
void ratcom_kevir_acc(double coef, double atomic_vir, double& val);
void ratcom_kevir_cu(double coef, double atomic_vir, double& val);


void rattle(time_prec dt, const pos_prec* xold, const pos_prec* yold,
            const pos_prec* zold);
void rattle_acc(time_prec, const pos_prec*, const pos_prec*, const pos_prec*);
void rattle_settle_acc(time_prec, const pos_prec*, const pos_prec*,
                       const pos_prec*);
void rattle_ch_acc(time_prec, const pos_prec*, const pos_prec*,
                   const pos_prec*);
// methylene and methyl groups
void rattle_methyl_cu(time_prec, const pos_prec*, const pos_prec*,
                      const pos_prec*);

void rattle2(time_prec dt, bool do_v);
void rattle2_acc(time_prec, bool);
void rattle2_settle_acc(time_prec, bool);
void rattle2_ch_acc(time_prec, bool);
void rattle2_methyl_cu(time_prec, bool);


void shake(time_prec dt, pos_prec* xnew, pos_prec* ynew, pos_prec* znew,
           const pos_prec* xold, const pos_prec* yold, const pos_prec* zold);
void shake_acc(time_prec dt, pos_prec* xnew, pos_prec* ynew, pos_prec* znew,
               const pos_prec* xold, const pos_prec* yold,
               const pos_prec* zold);
void shake_settle_acc(time_prec dt, pos_prec* xnew, pos_prec* ynew,
                      pos_prec* znew, const pos_prec* xold,
                      const pos_prec* yold, const pos_prec* zold);
void shake_ch_acc(time_prec, pos_prec*, pos_prec*, pos_prec*, const pos_prec*,
                  const pos_prec*, const pos_prec*);
void shake_methyl_cu(time_prec, pos_prec*, pos_prec*, pos_prec*,
                     const pos_prec*, const pos_prec*, const pos_prec*);
}
