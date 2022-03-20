#include "hippo/induce_donly.h"
#include "nblist.h"

namespace tinker {
void sparse_precond_build2() {}

void sparse_precond_apply2(const real (*rsd)[3], real (*zrsd)[3])
{
#if TINKER_CUDART
   if (ulist_version() & NBL_SPATIAL)
      sparse_precond_apply2_cu(rsd, zrsd);
   else
#endif
      sparse_precond_apply2_acc(rsd, zrsd);
}

void ulspred_save2(const real (*uind)[3])
{
   ulspred_save2_acc(uind);
}

void ulspred_sum2(real (*uind)[3])
{
   ulspred_sum2_acc(uind);
}
}
