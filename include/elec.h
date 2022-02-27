#pragma once
#include "mod.cflux.h"
#include "mod.charge.h"
#include "mod.chgpen.h"
#include "mod.chgpot.h"
#include "mod.ctrpot.h"
#include "mod.mpole.h"
#include "mod.polar.h"
#include "tool/rc_man.h"


namespace tinker {
bool use_ewald();


//====================================================================//


void pchg_data(rc_op);


//====================================================================//


// AMOBEA: multipole, polarization
// HIPPO: repulsion
void pole_data(rc_op);

void mdpuscale_data(rc_op);

void chgpen_data(rc_op op);
//====================================================================//


void elec_data(rc_op op);


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

bool amoebaplus_empole(int vers);
bool amoebaplus_epolar(int vers);
}
