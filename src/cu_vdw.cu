#include "e_vdw.h"
#include "md.h"

TINKER_NAMESPACE_BEGIN
extern void evdw_reduce_xyz();
extern void evdw_resolve_gradient();

void evdw_hal_cuda_impl_(int vers) {
  evdw_reduce_xyz();

  if (vers & calc::grad)
    evdw_resolve_gradient();
}
TINKER_NAMESPACE_END
