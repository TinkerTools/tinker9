#pragma once
#include "precision.h"

namespace tinker {
void diag_precond2(const real (*rsd)[3], real (*zrsd)[3]);

void sparse_precond_build2();
void sparse_precond_apply2(const real (*rsd)[3], real (*zrsd)[3]);
void sparse_precond_apply2_acc(const real (*)[3], real (*)[3]);
void sparse_precond_apply2_cu(const real (*)[3], real (*)[3]);

void ulspred_save2(const real (*uind)[3]);
void ulspred_sum2(real (*uind)[3]);
void ulspred_save2_acc(const real (*)[3]);
void ulspred_sum2_acc(real (*)[3]);
}
