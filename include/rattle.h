#pragma once
#include "mod.freeze.h"
#include "tool/rc_man.h"


namespace tinker {
bool use_rattle();
void rattle_data(rc_op);


void rattle(time_prec dt, const pos_prec* xold, const pos_prec* yold,
            const pos_prec* zold);
void rattle_acc(time_prec, const pos_prec*, const pos_prec*, const pos_prec*);
void rattle_settle_acc(time_prec, const pos_prec*, const pos_prec*,
                       const pos_prec*);

void rattle2(time_prec dt, bool do_v);
void rattle2_acc(time_prec, bool);
void rattle2_settle_acc(time_prec, bool);


void rattle_lp(time_prec dt, pos_prec* xnew, pos_prec* ynew, pos_prec* znew,
               const pos_prec* xold, const pos_prec* yold,
               const pos_prec* zold);
void rattle_acc_lp(time_prec dt, pos_prec* xnew, pos_prec* ynew, pos_prec* znew,
                   const pos_prec* xold, const pos_prec* yold,
                   const pos_prec* zold);
void rattle_settle_acc_lp(time_prec dt, pos_prec* xnew, pos_prec* ynew,
                          pos_prec* znew, const pos_prec* xold,
                          const pos_prec* yold, const pos_prec* zold);


void rattle2_lf(time_prec dt, pos_prec* vx_lp, pos_prec* vy_lp, pos_prec* vz_lp,
                const vel_prec* vx_new, const vel_prec* vy_new,
                const vel_prec* zx_new, const pos_prec* xold,
                const pos_prec* yold, const pos_prec* zold);
void rattle2_lf_acc(time_prec dt, pos_prec* vx_lp, pos_prec* vy_lp,
                    pos_prec* vz_lp, const vel_prec* vx_new,
                    const vel_prec* vy_new, const vel_prec* zx_new,
                    const pos_prec* xold, const pos_prec* yold,
                    const pos_prec* zold);
}
