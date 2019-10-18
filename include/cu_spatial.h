#pragma once
#include "spatial.h"
#include <thrust/execution_policy.h>
#include <thrust/extrema.h>
#include <thrust/sort.h>


TINKER_NAMESPACE_BEGIN
namespace cu_ker {
/**
 * \ingroup spatial
 * \brief A pair of `real` pointers on device.
 */
using RealPtrPair = thrust::pair<const real*, const real*>;


/**
 * \ingroup spatial
 * \brief Set the initial state of the K-D tree before space partitioning.
 * \param tree
 * Device pointer to the K-D tree.
 * \param x,y,z
 * Device pointers to the coordinates used to construct the neighbor list.
 * \param xm,ym,zm
 *
 */
__global__
void init_xyz_bbx_reorder(KDTree* RESTRICT tree, //
                          const real* RESTRICT x, const real* RESTRICT y,
                          const real* RESTRICT z, //
                          RealPtrPair xm, RealPtrPair ym, RealPtrPair zm)
{
   int n = tree->n_atom;
   int padded_n = tree->padded_n_atom;

   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
        i += blockDim.x * gridDim.x) {
      tree->xyz[i].x = x[i];
      tree->xyz[i].y = y[i];
      tree->xyz[i].z = z[i];
   }

   if (blockIdx.x == 0 && threadIdx.x == 0) {
      auto* RESTRICT bbx = tree->bbx;
      // bbx[0].begin = 0;
      bbx[0].end = tree->n_block;
      bbx[0].xmin = *xm.first;
      bbx[0].xmax = *xm.second;
      bbx[0].ymin = *ym.first;
      bbx[0].ymax = *ym.second;
      bbx[0].zmin = *zm.first;
      bbx[0].zmax = *zm.second;
   }

   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < padded_n;
        i += blockDim.x * gridDim.x) {
      tree->reorder[i] = i;
   }
}


__device__
inline int pick_axis(int l)
{
   return l - 3 * (l / 3);
}


__device__
inline int pow2(int l)
{
   return 1 << l;
}


__device__
inline int node_id(int l, int ml)
{
   return pow2(l) + ml - 1;
}


__device__
inline int left_child(int idx)
{
   return 2 * idx + 1;
}


struct x_less
{
   __device__
   bool operator()(const KDTree::Q& q1, const KDTree::Q& q2)
   {
      return q1.x < q2.x;
   }
};


struct y_less
{
   __device__
   bool operator()(const KDTree::Q& q1, const KDTree::Q& q2)
   {
      return q1.y < q2.y;
   }
};


struct z_less
{
   __device__
   bool operator()(const KDTree::Q& q1, const KDTree::Q& q2)
   {
      return q1.z < q2.z;
   }
};


__global__
void sort_level(int n, int l, int max_ml, KDTree* RESTRICT tree)
{
   auto* RESTRICT xyz = tree->xyz;
   auto* RESTRICT bbx = tree->bbx;
   auto* RESTRICT reorder = tree->reorder;
   auto* RESTRICT leafnode = tree->leaf_node;
   int ax = pick_axis(l);


   for (int ml = threadIdx.x + blockIdx.x * blockDim.x; ml < max_ml;
        ml += blockDim.x * gridDim.x) {
      int p = node_id(l, ml);
      int b0 = bbx[p].begin;
      int e0 = bbx[p].end;
      int b1 = b0 * KDTree::LEAF;
      int e1 = min(e0 * KDTree::LEAF, n);


      // if p does not have children
      if ((e0 - b0) == 1) {
         leafnode[b0] = p;
         thrust::sort_by_key(thrust::device, reorder + b1, reorder + e1,
                             xyz + b1);
         continue;
      } else if (e0 == b0) {
         continue;
      }


      // else p has children
      if (ax == 2)
         thrust::sort_by_key(thrust::device, xyz + b1, xyz + e1, reorder + b1,
                             z_less());
      else if (ax == 1)
         thrust::sort_by_key(thrust::device, xyz + b1, xyz + e1, reorder + b1,
                             y_less());
      else
         thrust::sort_by_key(thrust::device, xyz + b1, xyz + e1, reorder + b1,
                             x_less());


      int c0 = left_child(p);
      int c1 = c0 + 1;
      bbx[c0] = bbx[p];
      bbx[c1] = bbx[p];
      // Examples
      // [0, 10) -> 5; [0, 11) -> 5
      // [1, 10) -> 5; [1, 11) -> 6
      int cut = (b0 + e0) / 2;
      int cut1 = cut * KDTree::LEAF;
      int cut0 = cut1 - 1;
      // update bbx[c0] and bbx[c1]
      bbx[c0].end = cut;
      bbx[c1].begin = cut;
      if (ax == 2) {
         bbx[c0].zmax = xyz[cut0].z;
         bbx[c1].zmin = xyz[cut1].z;
      } else if (ax == 1) {
         bbx[c0].ymax = xyz[cut0].y;
         bbx[c1].ymin = xyz[cut1].y;
      } else {
         bbx[c0].xmax = xyz[cut0].x;
         bbx[c1].xmin = xyz[cut1].x;
      }
   }
}
}
TINKER_NAMESPACE_END
