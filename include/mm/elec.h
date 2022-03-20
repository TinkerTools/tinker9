#pragma once
#include "glob/cflux.h"
#include "glob/charge.h"
#include "glob/chgpen.h"
#include "glob/chgpot.h"
#include "glob/ctrpot.h"
#include "glob/mpole.h"
#include "glob/polar.h"
#include "tool/rcman.h"

namespace tinker {
bool use_ewald();

//====================================================================//

void pchg_data(RcOp);

//====================================================================//

// AMOBEA: multipole, polarization
// HIPPO: repulsion
void pole_data(RcOp);

void mdpuscale_data(RcOp);

void chgpen_data(RcOp);
//====================================================================//

void elec_data(RcOp);

//====================================================================//

void mpole_init(int vers);
void chkpole();
void rotpole();
void torque(int vers, grad_prec* dx, grad_prec* dy, grad_prec* dz);

void chkpole_acc();
void rotpole_acc();
void torque_acc(int vers, grad_prec* gx, grad_prec* gy, grad_prec* gz);
// void torque_cu(int vers);

//====================================================================//

bool amoeba_emplar(int vers);
bool amoeba_empole(int vers);
bool amoeba_epolar(int vers);

bool amoeba_echglj(int vers);
bool amoeba_echarge(int vers);
bool amoeba_evdw(int vers);

bool hippo_empole(int vers);
bool hippo_epolar(int vers);
}
