#include "cu_launch.h"
#include "cu_spatial.h"
#include "error.h"
#include "nblist.h"


TINKER_NAMESPACE_BEGIN
void spatial_data_init_cu(KDTreeUnit u, NBListUnit nu)
{
   int n = u->n_atom;
   int padded_n = u->padded_n_atom;
   int nblk = u->n_block;
   int nl = u->n_layer;

   auto* bbx = u->bbx;
   auto* reorder = u->reorder;
   auto* leafnode = u->leaf_node;

   int p2nl = pow2(nl);
   const real* x = nu->x;
   const real* y = nu->y;
   const real* z = nu->z;


   // find min/max of x/y/z
   auto xpair = thrust::minmax_element(thrust::device, x, x + n);
   auto ypair = thrust::minmax_element(thrust::device, y, y + n);
   auto zpair = thrust::minmax_element(thrust::device, z, z + n);
   // zero out bbx
   device_array::zero(p2nl, bbx);
   // initialize the tree
   launch_kernel1(padded_n, cu_ker::init_xyz_bbx_reorder, //
                  u.deviceptr(), x, y, z, xpair, ypair, zpair);


   // partition the space
   for (int l = 0; l < nl; ++l) {
      int max_ml = pow2(l);
      launch_kernel2(32, max_ml, cu_ker::sort_level, //
                     n, l, max_ml, u.deviceptr());
   }
}
TINKER_NAMESPACE_END
