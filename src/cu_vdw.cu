#include "e_vdw.h"
#include "md.h"
#include "nblist.h"

TINKER_NAMESPACE_BEGIN
template <int DO_G>
__device__ void //
ehal_pair_cuda(real rik, real& __restrict__ e, real& __restrict__ de) {}

void evdw_hal_cuda_impl_(int vers) {
  evdw_reduce_xyz();

  auto& vlst = *vlist_unit;
  real cut = vdw_switch_cut;
  real off = vdw_switch_off;
  real cut2 = cut * cut;
  real off2 = off * off;

  if (vers & calc::grad)
    evdw_resolve_gradient();
}
TINKER_NAMESPACE_END
