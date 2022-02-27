#pragma once
#include "macro.h"


namespace tinker {
void diag_precond3(const real (*rsd)[3], real (*zrsd)[3]);


void sparse_precond_build3();
void sparse_precond_apply3(const real (*rsd)[3], real (*zrsd)[3]);
void sparse_precond_apply3_acc(const real (*)[3], real (*)[3]);
void sparse_precond_apply3_cu(const real (*)[3], real (*)[3]);

void ulspred_save3(const real (*uind)[3]);
void ulspred_sum3(real (*uind)[3]);
void ulspred_save3_acc(const real (*)[3]);
void ulspred_sum3_acc(real (*)[3]);
}
