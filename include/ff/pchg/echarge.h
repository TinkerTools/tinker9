#pragma once
#include "mod/elecpchg.h"
#include "tool/rcman.h"

namespace tinker {

void echargeData(RcOp);

void echarge(int vers);

void echarge_nonewald(int);
void echarge_nonewald_acc(int);
void echarge_nonewald_cu(int);

void echarge_ewald_recip_self(int);
void echarge_ewald_fphi_self_acc(int);
void echarge_ewald_fphi_self_cu(int);

// void echarge_ewald_real(int);
void echarge_ewald_real_acc(int);
void echarge_ewald_real_cu(int);
}
