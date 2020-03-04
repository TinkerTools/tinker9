#pragma once
#include "macro.h"


TINKER_NAMESPACE_BEGIN
void diag_precond(const real (*rsd)[3], const real (*rsdp)[3], real (*zrsd)[3],
                  real (*zrsdp)[3]);


void sparse_precond_build();
void sparse_precond_apply(const real (*rsd)[3], const real (*rsdp)[3],
                          real (*zrsd)[3], real (*zrsdp)[3]);
void sparse_precond_apply_acc(const real (*)[3], const real (*)[3], real (*)[3],
                              real (*)[3]);
void sparse_precond_apply_cu(const real (*)[3], const real (*)[3], real (*)[3],
                             real (*)[3]);
TINKER_NAMESPACE_END
