#include "induce_aplus.h"
#include "nblist.h"


namespace tinker {
void sparse_precond_build3() {}


void sparse_precond_apply3(const real (*rsd)[3], real (*zrsd)[3])
{
#if TINKER_CUDART
   if (ulist_version() & NBL_SPATIAL)
      sparse_precond_apply3_cu(rsd, zrsd);
   else
#endif
      sparse_precond_apply3_acc(rsd, zrsd);
}


void ulspred_save3(const real (*uind)[3])
{
   ulspred_save3_acc(uind);
}


void ulspred_sum3(real (*uind)[3])
{
   ulspred_sum3_acc(uind);
}
}
