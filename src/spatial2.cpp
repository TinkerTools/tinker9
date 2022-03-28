#include "ff/spatial.h"
#include "math/maxmin.h"
#include "math/pow2.h"
#include "tool/darray.h"
#include "tool/error.h"

namespace tinker {
Spatial::~Spatial()
{
   darray::deallocate(iakpl, iak, lst);

   darray::deallocate(iakpl_rev, akpf, sorted, bnum);
   darray::deallocate(akc, half);
   darray::deallocate(update, xold, yold, zold);

   darray::deallocate(si1.bit0, si2.bit0, si3.bit0, si4.bit0);
}

void spatial2_cut(int& px, int& py, int& pz, int level);
void spatialDataAlloc(SpatialUnit& u, int n, double cutoff, double buffer, const real* x,
   const real* y, const real* z,                   //
   int nstype,                                     //
   int ns1, int (*js1)[2], int ns2, int (*js2)[2], //
   int ns3, int (*js3)[2], int ns4, int (*js4)[2])
{
   u = SpatialUnit::open();
   auto& st = *u;

   // output
   st.nakpl = 0;
   st.niak = 0;
   st.iakpl = nullptr;
   st.iak = nullptr;
   st.lst = nullptr;

   // internal
   st.n = n;
   st.nak = (n + Spatial::BLOCK - 1) / Spatial::BLOCK;
   st.nakp = (st.nak + 1) * st.nak / 2;
   st.nakpk = (st.nakp + Spatial::BLOCK - 1) / Spatial::BLOCK;
   int level = 1 + floorLog2(st.nak - 1);
   spatial2_cut(st.px, st.py, st.pz, level);
   st.cap_nakpl = 32 + 8 * st.nak;

   darray::allocate(st.cap_nakpl, &st.iakpl);
   darray::allocate(st.nak * Spatial::LSTCAP, &st.iak);
   darray::allocate(st.nak * Spatial::LSTCAP * 32, &st.lst);

   darray::allocate(st.nakp, &st.iakpl_rev);
   darray::allocate(st.nakpk, &st.akpf);
   darray::allocate(st.n, &st.sorted, &st.bnum);
   darray::allocate(st.nak, &st.akc, &st.half);

   darray::allocate(std::max(128, st.n * 2), &st.update);
   darray::allocate(st.n, &st.xold, &st.yold, &st.zold);

   st.fresh = 0;
   st.cutoff = cutoff;
   st.buffer = buffer;
   st.x = x;
   st.y = y;
   st.z = z;

   st.nstype = nstype;
   st.si1.init();
   st.si2.init();
   st.si3.init();
   st.si4.init();
   if (nstype >= 1) {
      st.si1.set(ns1, js1);
      darray::allocate(32 * st.cap_nakpl, &st.si1.bit0);
   }
   if (nstype >= 2) {
      st.si2.set(ns2, js2);
      darray::allocate(32 * st.cap_nakpl, &st.si2.bit0);
   }
   if (nstype >= 3) {
      st.si3.set(ns3, js3);
      darray::allocate(32 * st.cap_nakpl, &st.si3.bit0);
   }
   if (nstype >= 4) {
      st.si4.set(ns4, js4);
      darray::allocate(32 * st.cap_nakpl, &st.si4.bit0);
   }
}

void Spatial::ScaleInfo::init()
{
   js = nullptr;
   bit0 = nullptr;
   ns = 0;
}

void Spatial::ScaleInfo::set(int nns, int (*jjs)[2])
{
   ns = nns;
   js = jjs;
}

extern void spatial1_cut(int& px, int& py, int& pz, int level);
void spatial2_cut(int& px, int& py, int& pz, int level)
{
   spatial1_cut(px, py, pz, level);
   int pmax = maxOf(px, py, pz) + 1;
   pmax = std::min(pmax, 10);
   px = pmax;
   py = pmax;
   pz = pmax;
}
}
