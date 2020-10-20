#pragma once
#include "mod.cflux.h"
#include "mod.ctrpot.h"
#include "mod.mplpot.h"
#include "mod.mpole.h"
#include "tool/rc_man.h"


namespace tinker {
void cflux_data(rc_op op);
void alterchg();
void alterchg_acc();
void zero_pot();

void bndchg_acc1(real* pdelta);
void angchg_acc1(real* pdelta);
void dcflux_acc(int vers, grad_prec* gx, grad_prec* gy, grad_prec* gz,
            virial_buffer vir);

void dcflux(int vers, grad_prec* gx, grad_prec* gy, grad_prec* gz,
            virial_buffer vir);
}
