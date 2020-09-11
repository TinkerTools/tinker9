#pragma once
#include "macro.h"


namespace tinker {
void diag_precond2(const real (*rsd)[3], real (*zrsd)[3]);


void sparse_precond_build2();
void sparse_precond_apply2(const real (*rsd)[3], real (*zrsd)[3]);
void sparse_precond_apply_acc2(const real (*)[3], real (*)[3]);
void sparse_precond_apply_cu2(const real (*)[3], real (*)[3]);
}
