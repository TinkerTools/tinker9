#include "spatial2.h"
#include "tool/darray.h"
#include "tool/error.h"


namespace tinker {
Spatial2::~Spatial2()
{
   // output
   // internal
   darray::deallocate(sorted, bnum);
   darray::deallocate(update, xold, yold, zold);
}


void spatial2_cut(int& px, int& py, int& pz, int level)
{
   spatial_cut(px, py, pz, level);
   int pmax = max_of(px, py, pz);
   px = pmax + 1;
   py = pmax + 1;
   pz = pmax + 1;
}


void spatial2_data_alloc(Spatial2Unit& u, int n, double cutoff, double buffer,
                         const real* x, const real* y, const real* z, int nexcl,
                         int (*excl)[2], void* excl_scale, int NS)
{
   u = Spatial2Unit::open();
   auto& st = *u;


   // internal
   st.n = n;
   st.nak = (n + Spatial::BLOCK - 1) / Spatial::BLOCK;
   int level = 1 + floor_log2(st.nak - 1);
   spatial2_cut(st.px, st.py, st.pz, level);


   darray::allocate(st.n, &st.sorted, &st.bnum);
   darray::allocate(st.n * 2, &st.update);
   darray::allocate(st.n, &st.xold, &st.yold, &st.zold);


   st.rebuild = 1;
   st.cutoff = cutoff;
   st.buffer = buffer;
   st.x = x;
   st.y = y;
   st.z = z;


   u.update_deviceptr(st, WAIT_NEW_Q);
}


Spatial2Unit cspatial_v2_unit;
}
