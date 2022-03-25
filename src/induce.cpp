#include "ff/amoeba/induce.h"
#include "ff/nblist.h"

namespace tinker {
void sparse_precond_apply_acc(const real (*)[3], const real (*)[3], real (*)[3], real (*)[3]);
void sparse_precond_apply_cu(const real (*)[3], const real (*)[3], real (*)[3], real (*)[3]);

void ulspred_save_acc(const real (*)[3], const real (*)[3]);
void ulspred_sum_acc(real (*)[3], real (*)[3]);
}

namespace tinker {
void sparsePrecondBuild() {}

void sparsePrecondApply(
   const real (*rsd)[3], const real (*rsdp)[3], real (*zrsd)[3], real (*zrsdp)[3])
{
#if TINKER_CUDART
   if (ulist_version() & NBL_SPATIAL)
      sparse_precond_apply_cu(rsd, rsdp, zrsd, zrsdp);
   else
#endif
      sparse_precond_apply_acc(rsd, rsdp, zrsd, zrsdp);
}

void ulspredSave(const real (*uind)[3], const real (*uinp)[3])
{
   ulspred_save_acc(uind, uinp);
}

void ulspredSum(real (*uind)[3], real (*uinp)[3])
{
   ulspred_sum_acc(uind, uinp);
}
}
