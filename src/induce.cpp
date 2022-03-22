#include "ff/amoeba/induce.h"
#include "ff/nblist.h"

namespace tinker {
void sparse_precond_build() {}

void sparse_precond_apply(
   const real (*rsd)[3], const real (*rsdp)[3], real (*zrsd)[3], real (*zrsdp)[3])
{
#if TINKER_CUDART
   if (ulist_version() & NBL_SPATIAL)
      sparse_precond_apply_cu(rsd, rsdp, zrsd, zrsdp);
   else
#endif
      sparse_precond_apply_acc(rsd, rsdp, zrsd, zrsdp);
}

void ulspred_save(const real (*uind)[3], const real (*uinp)[3])
{
   ulspred_save_acc(uind, uinp);
}

void ulspred_sum(real (*uind)[3], real (*uinp)[3])
{
   ulspred_sum_acc(uind, uinp);
}
}
