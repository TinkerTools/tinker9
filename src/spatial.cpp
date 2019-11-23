#include "spatial.h"
#include "dev_array.h"
#include "error.h"
#include "md.h"
#include "nblist.h"


TINKER_NAMESPACE_BEGIN
Spatial::~Spatial()
{
   // output
   device_array::deallocate(lst, iak);
   // internal
   device_array::deallocate(sorted, boxnum);
   device_array::deallocate(naak, xakf, xakf_scan);
   device_array::deallocate(nearby);
   device_array::deallocate(ax_scan);
   device_array::deallocate(xkf);
   device_array::deallocate(update, xold, yold, zold);
}


void spatial_data_alloc(SpatialUnit& u, int n, double cutoff, double buffer,
                        const real* x, const real* y, const real* z)
{
   u = SpatialUnit::open();
   auto& st = *u;

   // output
   st.niak = 0;
   // internal
   st.n = n;
   st.nak = (n + Spatial::BLOCK - 1) / Spatial::BLOCK;
   int level = 1 + floor_log2(st.nak - 1);
   st.px = (level + 2) / 3;
   st.py = (level + 1) / 3;
   st.pz = (level + 0) / 3;
   st.nx = ct::pow2(st.px + st.py + st.pz);
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
   device_array::allocate(st.n, &st.update, &st.xold, &st.yold, &st.zold);

   st.rebuild = 1;
   st.cutoff = cutoff;
   st.buffer = buffer;
   st.x = x;
   st.y = y;
   st.z = z;

   u.update_deviceptr(st);
}


#if TINKER_HOST
void spatial_data_init_cu(SpatialUnit)
{
   TINKER_DUMMY_FUNCTION(spatial_data_init_cu);
}
#endif
TINKER_NAMESPACE_END
