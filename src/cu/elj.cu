#include "ff/evdw.h"
#include "ff/image.h"
#include "ff/spatial.h"
#include "ff/switch.h"
#include "math/switch.h"
#include "seq/add.h"
#include "seq/launch.h"
#include "seq/pair_lj.h"
#include "seq/triangle.h"

namespace tinker {
#include "elj_cu1.cc"

// special vdw14 interactions
template <class Ver>
__global__
void elj_cu2(CountBuffer restrict nebuf, EnergyBuffer restrict ebuf, VirialBuffer restrict vbuf, grad_prec* restrict gx,
   grad_prec* restrict gy, grad_prec* restrict gz, TINKER_IMAGE_PARAMS, real cut, real off, const real* restrict x,
   const real* restrict y, const real* restrict z, //
   int njvdw, const real* restrict radmin, const real* restrict epsilon, const int* restrict jvdw,
   const int* restrict mut, real vlam, Vdw vcouple, //
   real v4scale, int nvdw14, const int (*restrict vdw14ik)[2], const real* restrict radmin4,
   const real* restrict epsilon4)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;

   using ebuf_prec = EnergyBufferTraits::type;
   using vbuf_prec = VirialBufferTraits::type;
   ebuf_prec etl;
   vbuf_prec vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz;
   if CONSTEXPR (do_e) {
      etl = 0;
   }
   if CONSTEXPR (do_v) {
      vtlxx = 0;
      vtlyx = 0;
      vtlzx = 0;
      vtlyy = 0;
      vtlzy = 0;
      vtlzz = 0;
   }

   int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   for (int ii = ithread; ii < nvdw14; ii += blockDim.x * gridDim.x) {
      int i = vdw14ik[ii][0];
      int k = vdw14ik[ii][1];

      int jit = jvdw[i];
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      int imut = mut[i];

      int jkt = jvdw[k];
      real xr = xi - x[k];
      real yr = yi - y[k];
      real zr = zi - z[k];
      int kmut = mut[k];

      real r2 = image2(xr, yr, zr);
      if (r2 > off * off)
         continue;

      int pos = jit * njvdw + jkt;
      real rv = radmin[pos];
      real eps = epsilon[pos];
      real rv4 = radmin4[pos];
      real eps4 = epsilon4[pos];

      real r = REAL_SQRT(r2);
      real invr = REAL_RECIP(r);

      real e, de, e4, de4;
      real vlambda = pair_vlambda(vlam, vcouple, imut, kmut);
      pair_lj_v3<do_g, true, 0>(r, invr, vlambda, v4scale, rv, eps, cut, off, e, de);
      pair_lj_v3<do_g, true, 0>(r, invr, vlambda, v4scale, rv4, eps4, cut, off, e4, de4);
      e = e4 - e;
      if CONSTEXPR (do_g)
         de = de4 - de;

      if CONSTEXPR (do_e) {
         // if CONSTEXPR (do_a) {}
         etl += floatTo<ebuf_prec>(e);
      }
      if CONSTEXPR (do_g) {
         de *= invr;
         real dedx = de * xr;
         real dedy = de * yr;
         real dedz = de * zr;
         if CONSTEXPR (do_v) {
            vtlxx += floatTo<vbuf_prec>(xr * dedx);
            vtlyx += floatTo<vbuf_prec>(yr * dedx);
            vtlzx += floatTo<vbuf_prec>(zr * dedx);
            vtlyy += floatTo<vbuf_prec>(yr * dedy);
            vtlzy += floatTo<vbuf_prec>(zr * dedy);
            vtlzz += floatTo<vbuf_prec>(zr * dedz);
         }
         atomic_add(dedx, gx, i);
         atomic_add(dedy, gy, i);
         atomic_add(dedz, gz, i);
         atomic_add(-dedx, gx, k);
         atomic_add(-dedy, gy, k);
         atomic_add(-dedz, gz, k);
      }
   }

   if CONSTEXPR (do_e) {
      atomic_add(etl, ebuf, ithread);
   }
   if CONSTEXPR (do_v) {
      atomic_add(vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz, vbuf, ithread);
   }
}

template <class Ver>
static void elj_cu4()
{
   const auto& st = *cspatial_v2_unit;
   const real cut = switchCut(Switch::VDW);
   const real off = switchOff(Switch::VDW);

   int ngrid = gpuGridSize(BLOCK_DIM);
   elj_cu1<Ver><<<ngrid, BLOCK_DIM, 0, g::s0>>>(st.n, TINKER_IMAGE_ARGS, nev, ev, vir_ev, devx, devy, devz, cut, off,
      st.si2.bit0, nvexclude, vexclude, vexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl, st.niak,
      st.iak, st.lst, njvdw, radmin, epsilon, jvdw, mut, vlam, vcouple);
}

template <class Ver>
static void elj_cu5()
{
   if (nvdw14 > 0) {
      const auto& st = *cspatial_v2_unit;
      const real cut = switchCut(Switch::VDW);
      const real off = switchOff(Switch::VDW);

      launch_k1b(g::s0, nvdw14, elj_cu2<Ver>,                                              //
         nev, ev, vir_ev, devx, devy, devz, TINKER_IMAGE_ARGS, cut, off, st.x, st.y, st.z, //
         njvdw, radmin, epsilon, jvdw, mut, vlam, vcouple,                                 //
         v4scale, nvdw14, vdw14ik, radmin4, epsilon4);
   }
}

void elj14_cu(int vers)
{
   if (vers == calc::v0)
      elj_cu5<calc::V0>();
   else if (vers == calc::v1)
      elj_cu5<calc::V1>();
   else if (vers == calc::v3)
      elj_cu5<calc::V3>();
   else if (vers == calc::v4)
      elj_cu5<calc::V4>();
   else if (vers == calc::v5)
      elj_cu5<calc::V5>();
   else if (vers == calc::v6)
      elj_cu5<calc::V6>();
}

void elj_cu(int vers)
{
   if (vers == calc::v0)
      elj_cu4<calc::V0>();
   else if (vers == calc::v1)
      elj_cu4<calc::V1>();
   else if (vers == calc::v3)
      elj_cu4<calc::V3>();
   else if (vers == calc::v4)
      elj_cu4<calc::V4>();
   else if (vers == calc::v5)
      elj_cu4<calc::V5>();
   else if (vers == calc::v6)
      elj_cu4<calc::V6>();

   elj14_cu(vers);
}
}
