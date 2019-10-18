#include "spatial.h"
#include "dev_array.h"
#include "error.h"
#include "md.h"
#include "nblist.h"


TINKER_NAMESPACE_BEGIN
KDTree::~KDTree()
{
   device_array::deallocate(xyz, bbx, reorder, leaf_node);
}


static void spatial_data_alloc(KDTreeUnit& u, int n)
{
   u = KDTreeUnit::open();
   auto& st = *u;

   st.n_atom = n;
   st.n_block = KDTree::nleaf(st.n_atom);
   st.padded_n_atom = st.n_block * KDTree::LEAF;
   st.n_layer = KDTree::nlayer(st.n_atom);

   device_array::allocate(st.n_atom, &st.xyz);
   device_array::allocate(pow2(st.n_layer), &st.bbx);
   device_array::allocate(st.padded_n_atom, &st.reorder);
   device_array::allocate(st.n_block, &st.leaf_node);

   u.update_deviceptr(st);
}


extern void spatial_data_init_cu(KDTreeUnit, NBListUnit);
void spatial_data(rc_op op)
{
   if (op & rc_dealloc)
      KDTreeUnit::clear();

   if (op & rc_alloc) {
      if (vlist_unit.valid())
         spatial_data_alloc(vtree_unit, n);

      if (mlist_unit.valid())
         spatial_data_alloc(mtree_unit, n);
   }

   if (op & rc_init) {
      if (vlist_unit.valid())
         spatial_data_init_cu(vtree_unit, vlist_unit);

      if (mlist_unit.valid())
         spatial_data_init_cu(mtree_unit, mlist_unit);
   }
}


#if !TINKER_CUDART
void spatial_data_init_cu(KDTreeUnit, NBListUnit)
{
   TINKER_THROW("This dummy function should not have been called.");
}
#endif
TINKER_NAMESPACE_END
