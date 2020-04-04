#include "induce.h"
#include "nblist.h"


TINKER_NAMESPACE_BEGIN
void sparse_precond_build() {}


void sparse_precond_apply(const real (*rsd)[3], const real (*rsdp)[3],
                          real (*zrsd)[3], real (*zrsdp)[3])
{
#if TINKER_CUDART
   if (ulist_version() & NBL_SPATIAL)
      sparse_precond_apply_cu(rsd, rsdp, zrsd, zrsdp);
   else
#endif
      sparse_precond_apply_acc(rsd, rsdp, zrsd, zrsdp);
}
TINKER_NAMESPACE_END
