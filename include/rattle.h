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


void shake(time_prec dt, pos_prec* xnew, pos_prec* ynew, pos_prec* znew,
           const pos_prec* xold, const pos_prec* yold, const pos_prec* zold);
void shake_acc(time_prec dt, pos_prec* xnew, pos_prec* ynew, pos_prec* znew,
               const pos_prec* xold, const pos_prec* yold,
               const pos_prec* zold);
void shake_settle_acc(time_prec dt, pos_prec* xnew, pos_prec* ynew,
                      pos_prec* znew, const pos_prec* xold,
                      const pos_prec* yold, const pos_prec* zold);

void shake2(time_prec dt, const vel_prec* vxold, const vel_prec* vyold,
            const vel_prec* vzold, const vel_prec* vxnew, const vel_prec* vynew,
            const vel_prec* vznew, const pos_prec* xold, const pos_prec* yold,
            const pos_prec* zold);
void shake2_acc(time_prec dt, const vel_prec* vxold, const vel_prec* vyold,
                const vel_prec* vzold, const vel_prec* vxnew,
                const vel_prec* vynew, const vel_prec* vznew,
                const pos_prec* xold, const pos_prec* yold,
                const pos_prec* zold);
}
