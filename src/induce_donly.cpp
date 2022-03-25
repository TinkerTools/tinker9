#include "ff/hippo/induce.h"
#include "ff/nblist.h"

namespace tinker {
void sparse_precond_apply2_acc(const real (*)[3], real (*)[3]);
void sparse_precond_apply2_cu(const real (*)[3], real (*)[3]);
void sparsePrecondBuild2() {}

void sparsePrecondApply2(const real (*rsd)[3], real (*zrsd)[3])
{
#if TINKER_CUDART
   if (ulist_version() & NBL_SPATIAL)
      sparse_precond_apply2_cu(rsd, zrsd);
   else
#endif
      sparse_precond_apply2_acc(rsd, zrsd);
}

void ulspred_save2_acc(const real (*)[3]);
void ulspred_sum2_acc(real (*)[3]);
void ulspredSave2(const real (*uind)[3])
{
   ulspred_save2_acc(uind);
}

void ulspredSum2(real (*uind)[3])
{
   ulspred_sum2_acc(uind);
}
}
