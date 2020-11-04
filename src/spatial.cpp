#include "spatial.h"
#include "md.h"
#include "nblist.h"
#include "tool/darray.h"
#include "tool/error.h"


namespace tinker {
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


void spatial_cut_v1(int& px, int& py, int& pz, int level)
{
   px = (level + 2) / 3;
   py = (level + 1) / 3;
   pz = (level + 0) / 3;
}
void spatial_cut_v2(int& px, int& py, int& pz, int level)
{
   px = (level + 0) / 3;
   py = (level + 1) / 3;
   pz = (level + 2) / 3;
}
void spatial_cut_v3(int& px, int& py, int& pz, int level)
{
   // triclinic frac(1,1,1) -> cart(x,y,z)
   // x = (fz * l1.z + fy * l1.y + fx * l1.x)
   // y = (fz * l2.z + fy * l2.y)
   // z = (fz * l3.z)


   double3 l1 = make_double3(lvec1.x, lvec1.y, lvec1.z);
   double3 l2 = make_double3(lvec2.x, lvec2.y, lvec2.z);
   double3 l3 = make_double3(lvec3.x, lvec3.y, lvec3.z);
   px = 0;
   py = 0;
   pz = 0;
   const double ratio = 0.95;
   for (int i = 0; i < level; ++i) {
      double xx = l1.z + l1.y + l1.x;
      double yy = l2.z + l2.y;
      double zz = l3.z;


      if ((zz > ratio * xx) && (zz > ratio * yy)) {
         // if z is approximately the longest, cut c-axis by half
         l1.z /= 2;
         l2.z /= 2;
         l3.z /= 2;
         pz += 1;
      } else if (yy > ratio * xx) {
         // if y is approximately the longest, cut b-axis by half
         l1.y /= 2;
         l2.y /= 2;
         l3.y /= 2;
         py += 1;
      } else {
         // if x is longest, cut a-axis by half
         l1.x /= 2;
         l2.x /= 2;
         l3.x /= 2;
         px += 1;
      }
   }
}
void spatial_cut(int& px, int& py, int& pz, int level)
{
   spatial_cut_v3(px, py, pz, level);
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
   spatial_cut(st.px, st.py, st.pz, level);
   st.nx = pow2(st.px + st.py + st.pz);
   st.nxk = (st.nx + Spatial::BLOCK - 1) / Spatial::BLOCK;
   st.near = 0;
   st.xak_sum = 0;
   st.iak_cap = 0;

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

   st.fresh = 0;
   st.cutoff = cutoff;
   st.buffer = buffer;
   st.x = x;
   st.y = y;
   st.z = z;

   u.update_deviceptr(st, g::q0);
   wait_for(g::q0);
}
}
