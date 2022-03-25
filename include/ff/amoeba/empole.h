#pragma once
#include "ff/elec.h"
#include "mod/elecamoeba.h"

namespace tinker {
void mpoleInit(int vers);
void torque(int vers, grad_prec* dx, grad_prec* dy, grad_prec* dz);

void empoleData(RcOp);
void empole(int vers);

void empole_nonewald(int vers);
void empole_ewald(int vers);
void empole_ewald_real_self(int vers);
void empole_ewald_recip(int vers);

void empole_nonewald_acc(int vers);
void empole_ewald_recip_acc(int vers);
void empole_ewald_real_self_acc(int vers);
void empole_nonewald_cu(int vers);
void empole_ewald_real_self_cu(int vers);

template <class Ver, int CFLX>
void empole_generic_ewald_recip_acc();
}
