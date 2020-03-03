#include "spatial.h"
#include "darray.h"
#include "error.h"
#include "md.h"
#include "nblist.h"


TINKER_NAMESPACE_BEGIN
Spatial::~Spatial()
{
   // output
   darray::deallocate(lst, iak);
   // internal
   darray::deallocate(sorted, boxnum);
   darray::deallocate(naak, xakf, xakf_scan);
   darray::deallocate(nearby);
   darray::deallocate(ax_scan);
   darray::deallocate(xkf);
   darray::deallocate(update, xold, yold, zold);
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
   st.px = (level + 0) / 3;
   st.py = (level + 1) / 3;
   st.pz = (level + 2) / 3;
   st.nx = pow2(st.px + st.py + st.pz);
   st.nxk = (st.nx + Spatial::BLOCK - 1) / Spatial::BLOCK;
   st.near = 0;
   st.xak_sum = 0;
   st.xak_sum_cap = 0;

   // output
   st.iak = nullptr;
   st.lst = nullptr;
   // internal
   darray::allocate(st.n, &st.sorted, &st.boxnum);
   darray::allocate(st.nak, &st.naak, &st.xakf, &st.xakf_scan);
   darray::allocate(st.nx, &st.nearby);
   darray::allocate(st.nx + 1, &st.ax_scan);
   darray::allocate(st.nak * st.nxk, &st.xkf);
   darray::allocate(st.n, &st.update, &st.xold, &st.yold, &st.zold);

   st.rebuild = 1;
   st.cutoff = cutoff;
   st.buffer = buffer;
   st.x = x;
   st.y = y;
   st.z = z;

   u.update_deviceptr(st, WAIT_NEW_Q);
}
TINKER_NAMESPACE_END
