#pragma once
#include "macro.h"


namespace tinker {
void diag_precond(const real (*rsd)[3], const real (*rsdp)[3], real (*zrsd)[3],
                  real (*zrsdp)[3]);


void sparse_precond_build();
void sparse_precond_apply(const real (*rsd)[3], const real (*rsdp)[3],
                          real (*zrsd)[3], real (*zrsdp)[3]);
void sparse_precond_apply_acc(const real (*)[3], const real (*)[3], real (*)[3],
                              real (*)[3]);
void sparse_precond_apply_cu(const real (*)[3], const real (*)[3], real (*)[3],
                             real (*)[3]);


void ulspred_save(const real (*uind)[3], const real (*uinp)[3]);
void ulspred_sum(real (*uind)[3], real (*uinp)[3]);
void ulspred_save_acc(const real (*)[3], const real (*)[3]);
void ulspred_sum_acc(real (*)[3], real (*)[3]);
}
