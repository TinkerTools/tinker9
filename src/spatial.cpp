#include "spatial.h"
#include "dev_array.h"
#include "error.h"
#include "md.h"
#include "nblist.h"


TINKER_NAMESPACE_BEGIN
Spatial::~Spatial()
{
   // output
   device_array::deallocate(lst); // iak and lst were allocated together
   // internal
   device_array::deallocate(sorted, boxnum);
   device_array::deallocate(naak, xakf, xakf_scan);
   device_array::deallocate(nearby);
   device_array::deallocate(ax_scan);
   device_array::deallocate(xkf);
}


static void spatial_data_alloc(SpatialUnit& u, int n)
{
   u = SpatialUnit::open();
   auto& st = *u;

   // output
   st.niak = 0;
   // internal
   st.n = n;
   st.nak = (n + Spatial::BLOCK - 1) / Spatial::BLOCK;
   int level = 1 + builtin_floor_log2(st.nak - 1);
   st.px = (level + 2) / 3;
   st.py = (level + 1) / 3;
   st.pz = (level + 0) / 3;
   st.nx = pow2(st.px + st.py + st.pz);
   st.nxk = (st.nx + Spatial::BLOCK - 1) / Spatial::BLOCK;
   st.near = 0;
   st.xak_sum = 0;
   st.xak_sum_cap = 0;

   // output
   st.iak = nullptr;
   st.lst = nullptr;
   // internal
   device_array::allocate(st.n, &st.sorted, &st.boxnum);
   device_array::allocate(st.nak, &st.naak, &st.xakf, &st.xakf_scan);
   device_array::allocate(st.nx, &st.nearby);
   device_array::allocate(st.nx + 1, &st.ax_scan);
   device_array::allocate(st.nak * st.nxk, &st.xkf);

   u.update_deviceptr(st);
}


extern void spatial_data_init_cu(SpatialUnit, NBListUnit);
extern void thrust_cache_dealloc();
extern void thrust_cache_alloc();
void spatial_data(rc_op op)
{
   if (op & rc_dealloc) {
      SpatialUnit::clear();
      thrust_cache_dealloc();
      vspatial_unit.close();
      mspatial_unit.close();
   }


   if (op & rc_alloc) {
      bool alloc_thrust_cache = false;

      if (vlist_unit.valid()) {
         spatial_data_alloc(vspatial_unit, n);
         alloc_thrust_cache = true;
      }

      if (mlist_unit.valid()) {
         spatial_data_alloc(mspatial_unit, n);
         alloc_thrust_cache = true;
      }

      if (alloc_thrust_cache)
         thrust_cache_alloc();
   }


   if (op & (rc_init | rc_evolve)) {
      if (vlist_unit.valid()) {
         // evdw_reduce_xyz(); // if is not always called after nblist_data()
         spatial_data_init_cu(vspatial_unit, vlist_unit);
      }

      if (mlist_unit.valid())
         spatial_data_init_cu(mspatial_unit, mlist_unit);
   }
}


#if TINKER_HOST
void spatial_data_init_cu(SpatialUnit, NBListUnit)
{
   TINKER_THROW("This dummy function should not have been called.");
}


void thrust_cache_dealloc()
{
   spatial_data_init_cu(0, 0);
}


void thrust_cache_alloc()
{
   spatial_data_init_cu(0, 0);
}
#endif
TINKER_NAMESPACE_END
