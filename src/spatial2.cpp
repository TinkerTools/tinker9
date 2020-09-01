#include "spatial2.h"
#include "tool/darray.h"
#include "tool/error.h"
#include <algorithm>
#include <set>


namespace tinker {
Spatial2::~Spatial2()
{
   darray::deallocate(iakpl, iak, lst);


   darray::deallocate(iakpl_rev, akpf, iakbuf, sorted, bnum);
   darray::deallocate(akc, half);
   darray::deallocate(update, xold, yold, zold);


   darray::deallocate(si1.bit0, si1.bit1, si2.bit0, si2.bit1, //
                      si3.bit0, si3.bit1, si4.bit0, si4.bit1);
}


void spatial2_data_alloc(
   Spatial2Unit& u, int n, double cutoff, double buffer, const real* x,
   const real* y, const real* z,                                      //
   int nstype,                                                        //
   int ns1, int (*js1)[2], real* ks1, const std::vector<double>& vs1, //
   int ns2, int (*js2)[2], real* ks2, const std::vector<double>& vs2, //
   int ns3, int (*js3)[2], real* ks3, const std::vector<double>& vs3, //
   int ns4, int (*js4)[2], real* ks4, const std::vector<double>& vs4)
{
   u = Spatial2Unit::open();
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
   int level = 1 + floor_log2(st.nak - 1);
   spatial2_cut(st.px, st.py, st.pz, level);
   st.cap_nakpl = 8 * st.nak;


   darray::allocate(st.cap_nakpl, &st.iakpl);
   darray::allocate(st.nak * Spatial2::LSTCAP, &st.iak);
   darray::allocate(st.nak * Spatial2::LSTCAP * 32, &st.lst);


   darray::allocate(st.nakp, &st.iakpl_rev);
   darray::allocate(st.nakpk, &st.akpf);
   darray::allocate(st.nak * Spatial2::LSTCAP, &st.iakbuf);
   darray::allocate(st.n, &st.sorted, &st.bnum);
   darray::allocate(st.nak, &st.akc, &st.half);


   darray::allocate(std::max(128, st.n * 2), &st.update);
   darray::allocate(st.n, &st.xold, &st.yold, &st.zold);


   st.rebuild = 1;
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
      st.si1.set(ns1, js1, ks1, vs1);
      darray::allocate(32 * st.cap_nakpl, &st.si1.bit0, &st.si1.bit1);
   }
   if (nstype >= 2) {
      st.si2.set(ns2, js2, ks2, vs2);
      darray::allocate(32 * st.cap_nakpl, &st.si2.bit0, &st.si2.bit1);
   }
   if (nstype >= 3) {
      st.si3.set(ns3, js3, ks3, vs3);
      darray::allocate(32 * st.cap_nakpl, &st.si3.bit0, &st.si3.bit1);
   }
   if (nstype >= 4) {
      st.si4.set(ns4, js4, ks4, vs4);
      darray::allocate(32 * st.cap_nakpl, &st.si4.bit0, &st.si4.bit1);
   }
}


void Spatial2::ScaleInfo::init()
{
   js = nullptr;
   ks = nullptr;
   bit0 = nullptr;
   bit1 = nullptr;
   sc0 = 0;
   sc1 = 0;
   sc2 = 0;
   sc3 = 0;
   ns = 0;
}


void Spatial2::ScaleInfo::set(int nns, int (*jjs)[2], real* kks,
                              const std::vector<double>& vs)
{
   ns = nns;
   js = jjs;
   ks = kks;


   std::vector<int> vsint;
   vsint.push_back(1 * MULT);
   for (auto fs : vs)
      vsint.push_back(fs * MULT);
   std::set<int> ss(vsint.begin(), vsint.end());
   if (ss.size() > CAPACITY) {
      TINKER_THROW(format(
         "Spatial2::ScaleInfo::set -- more than %zu 1-x scale factors are used",
         ss.size()));
   }


   std::set<real> sfloat;
   sfloat.insert(1);
   for (auto fs : vs)
      sfloat.insert(fs);
   if (sfloat.size() != ss.size()) {
      TINKER_THROW(
         "Spatial2::ScaleInfo::set -- error in parsing the scale factors");
   }


   std::vector<real> v2(sfloat.begin(), sfloat.end());
   std::reverse(v2.begin(), v2.end());
   if (v2.size() > 0)
      sc0 = vs[0];
   if (v2.size() > 1)
      sc1 = vs[1];
   if (v2.size() > 2)
      sc2 = vs[2];
   if (v2.size() > 3)
      sc3 = vs[3];
}


void spatial2_cut(int& px, int& py, int& pz, int level)
{
   spatial_cut(px, py, pz, level);
   int pmax = max_of(px, py, pz) + 1;
   pmax = std::min(pmax, 10);
   px = pmax;
   py = pmax;
   pz = pmax;
}


Spatial2Unit cspatial_v2_unit;
}
