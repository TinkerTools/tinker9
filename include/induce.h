#pragma once
#include "e_polar.h"
#include "md.h"
#include "nblist.h"
#include "seq_damp.h"
#include "seq_image.h"


TINKER_NAMESPACE_BEGIN
void diag_precond(const real (*rsd)[3], const real (*rsdp)[3], real (*zrsd)[3],
                  real (*zrsdp)[3]);


void sparse_precond_build();
void sparse_precond_apply(const real (*rsd)[3], const real (*rsdp)[3],
                          real (*zrsd)[3], real (*zrsdp)[3]);
TINKER_NAMESPACE_END
