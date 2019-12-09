#include "induce.h"


TINKER_NAMESPACE_BEGIN
void sparse_precond_build()
{
   extern void sparse_precond_build_acc();
#if TINKER_CUDART
   if (ulist_version() == NBList::spatial) {
      // do nothing
   } else
#endif
      sparse_precond_build_acc();
}


void sparse_precond_apply(const real (*rsd)[3], const real (*rsdp)[3],
                          real (*zrsd)[3], real (*zrsdp)[3])
{
   extern void sparse_precond_apply_acc(const real(*)[3], const real(*)[3],
                                        real(*)[3], real(*)[3]);
#if TINKER_CUDART
   if (ulist_version() == NBList::spatial) {
      extern void sparse_precond_apply_cu(const real(*)[3], const real(*)[3],
                                          real(*)[3], real(*)[3]);
      sparse_precond_apply_cu(rsd, rsdp, zrsd, zrsdp);
   } else
#endif
      sparse_precond_apply_acc(rsd, rsdp, zrsd, zrsdp);
}
TINKER_NAMESPACE_END
