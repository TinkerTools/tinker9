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
// ck.py Version 2.0.3
template <class Ver>
__global__
void elj_cu1(int n, TINKER_IMAGE_PARAMS, CountBuffer restrict nev, EnergyBuffer restrict ev,
   VirialBuffer restrict vev, grad_prec* restrict gx, grad_prec* restrict gy,
   grad_prec* restrict gz, real cut, real off, const unsigned* restrict info, int nexclude,
   const int (*restrict exclude)[2], const real* restrict exclude_scale, const real* restrict x,
   const real* restrict y, const real* restrict z, const Spatial::SortedAtom* restrict sorted,
   int nakpl, const int* restrict iakpl, int niak, const int* restrict iak, const int* restrict lst,
   int njvdw, const real* restrict radmin, const real* restrict epsilon, const int* restrict jvdw,
   const int* restrict mut, real vlam, Vdw vcouple)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);

   int nevtl;
   if CONSTEXPR (do_a) {
      nevtl = 0;
   }
   using ebuf_prec = EnergyBufferTraits::type;
   ebuf_prec evtl;
   if CONSTEXPR (do_e) {
      evtl = 0;
   }
   using vbuf_prec = VirialBufferTraits::type;
   vbuf_prec vevtlxx, vevtlyx, vevtlzx, vevtlyy, vevtlzy, vevtlzz;
   if CONSTEXPR (do_v) {
      vevtlxx = 0;
      vevtlyx = 0;
      vevtlzx = 0;
      vevtlyy = 0;
      vevtlzy = 0;
      vevtlzz = 0;
   }
   real xi;
   real yi;
   real zi;
   real xk;
   real yk;
   real zk;
   real fix;
   real fiy;
   real fiz;
   real fkx;
   real fky;
   real fkz;
   int ijvdw;
   int imut;
   int kjvdw;
   int kmut;

   //* /
   for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
      if CONSTEXPR (do_g) {
         fix = 0;
         fiy = 0;
         fiz = 0;
         fkx = 0;
         fky = 0;
         fkz = 0;
      }

      int i = exclude[ii][0];
      int k = exclude[ii][1];
      real scalea = exclude_scale[ii];

      xi = x[i];
      yi = y[i];
      zi = z[i];
      xk = x[k];
      yk = y[k];
      zk = z[k];
      ijvdw = jvdw[i];
      imut = mut[i];
      kjvdw = jvdw[k];
      kmut = mut[k];

      constexpr bool incl = true;
      real xr = xi - xk;
      real yr = yi - yk;
      real zr = zi - zk;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         real r = REAL_SQRT(r2);
         real invr = REAL_RECIP(r);
         real rv = radmin[ijvdw * njvdw + kjvdw];
         real eps = epsilon[ijvdw * njvdw + kjvdw];
         real e, de;
         real vlambda = pair_vlambda(vlam, vcouple, imut, kmut);
         pair_lj_v3<do_g, true, 0>(r, invr, vlambda, scalea, rv, eps, cut, off, e, de);
         if CONSTEXPR (do_e) {
            evtl += floatTo<ebuf_prec>(e);
            if CONSTEXPR (do_a) {
               if (scalea != 0 and e != 0)
                  nevtl += 1;
            }
         }
         if CONSTEXPR (do_g) {
            real dedx, dedy, dedz;
            de = de * invr;
            dedx = de * xr;
            dedy = de * yr;
            dedz = de * zr;
            fix += dedx;
            fiy += dedy;
            fiz += dedz;
            fkx -= dedx;
            fky -= dedy;
            fkz -= dedz;
            if CONSTEXPR (do_v) {
               vevtlxx += floatTo<vbuf_prec>(xr * dedx);
               vevtlyx += floatTo<vbuf_prec>(yr * dedx);
               vevtlzx += floatTo<vbuf_prec>(zr * dedx);
               vevtlyy += floatTo<vbuf_prec>(yr * dedy);
               vevtlzy += floatTo<vbuf_prec>(zr * dedy);
               vevtlzz += floatTo<vbuf_prec>(zr * dedz);
            }
         }
      }

      if CONSTEXPR (do_g) {
         atomic_add(fix, gx, i);
         atomic_add(fiy, gy, i);
         atomic_add(fiz, gz, i);
         atomic_add(fkx, gx, k);
         atomic_add(fky, gy, k);
         atomic_add(fkz, gz, k);
      }
   }
   // */

   for (int iw = iwarp; iw < nakpl; iw += nwarp) {
      if CONSTEXPR (do_g) {
         fix = 0;
         fiy = 0;
         fiz = 0;
         fkx = 0;
         fky = 0;
         fkz = 0;
      }

      int tri, tx, ty;
      tri = iakpl[iw];
      tri_to_xy(tri, tx, ty);

      int iid = ty * WARP_SIZE + ilane;
      int atomi = min(iid, n - 1);
      int i = sorted[atomi].unsorted;
      int kid = tx * WARP_SIZE + ilane;
      int atomk = min(kid, n - 1);
      int k = sorted[atomk].unsorted;
      xi = sorted[atomi].x;
      yi = sorted[atomi].y;
      zi = sorted[atomi].z;
      xk = sorted[atomk].x;
      yk = sorted[atomk].y;
      zk = sorted[atomk].z;

      ijvdw = jvdw[i];
      imut = mut[i];
      kjvdw = jvdw[k];
      kmut = mut[k];

      unsigned int info0 = info[iw * WARP_SIZE + ilane];
      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         bool incl = iid < kid and kid < n;
         int srcmask = 1 << srclane;
         incl = incl and (info0 & srcmask) == 0;
         real xr = xi - xk;
         real yr = yi - yk;
         real zr = zi - zk;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real r = REAL_SQRT(r2);
            real invr = REAL_RECIP(r);
            real rv = radmin[ijvdw * njvdw + kjvdw];
            real eps = epsilon[ijvdw * njvdw + kjvdw];
            real e, de;
            real vlambda = pair_vlambda(vlam, vcouple, imut, kmut);
            pair_lj_v3<do_g, true, 1>(r, invr, vlambda, 1, rv, eps, cut, off, e, de);
            if CONSTEXPR (do_e) {
               evtl += floatTo<ebuf_prec>(e);
               if CONSTEXPR (do_a) {
                  if (e != 0)
                     nevtl += 1;
               }
            }
            if CONSTEXPR (do_g) {
               real dedx, dedy, dedz;
               de = de * invr;
               dedx = de * xr;
               dedy = de * yr;
               dedz = de * zr;
               fix += dedx;
               fiy += dedy;
               fiz += dedz;
               fkx -= dedx;
               fky -= dedy;
               fkz -= dedz;
               if CONSTEXPR (do_v) {
                  vevtlxx += floatTo<vbuf_prec>(xr * dedx);
                  vevtlyx += floatTo<vbuf_prec>(yr * dedx);
                  vevtlzx += floatTo<vbuf_prec>(zr * dedx);
                  vevtlyy += floatTo<vbuf_prec>(yr * dedy);
                  vevtlzy += floatTo<vbuf_prec>(zr * dedy);
                  vevtlzz += floatTo<vbuf_prec>(zr * dedz);
               }
            }
         }

         iid = __shfl_sync(ALL_LANES, iid, ilane + 1);
         xi = __shfl_sync(ALL_LANES, xi, ilane + 1);
         yi = __shfl_sync(ALL_LANES, yi, ilane + 1);
         zi = __shfl_sync(ALL_LANES, zi, ilane + 1);
         ijvdw = __shfl_sync(ALL_LANES, ijvdw, ilane + 1);
         imut = __shfl_sync(ALL_LANES, imut, ilane + 1);
         if CONSTEXPR (do_g) {
            fix = __shfl_sync(ALL_LANES, fix, ilane + 1);
            fiy = __shfl_sync(ALL_LANES, fiy, ilane + 1);
            fiz = __shfl_sync(ALL_LANES, fiz, ilane + 1);
         }
      }

      if CONSTEXPR (do_g) {
         atomic_add(fix, gx, i);
         atomic_add(fiy, gy, i);
         atomic_add(fiz, gz, i);
         atomic_add(fkx, gx, k);
         atomic_add(fky, gy, k);
         atomic_add(fkz, gz, k);
      }
   }

   for (int iw = iwarp; iw < niak; iw += nwarp) {
      if CONSTEXPR (do_g) {
         fix = 0;
         fiy = 0;
         fiz = 0;
         fkx = 0;
         fky = 0;
         fkz = 0;
      }

      int ty = iak[iw];
      int atomi = ty * WARP_SIZE + ilane;
      int i = sorted[atomi].unsorted;
      int atomk = lst[iw * WARP_SIZE + ilane];
      int k = sorted[atomk].unsorted;
      xi = sorted[atomi].x;
      yi = sorted[atomi].y;
      zi = sorted[atomi].z;
      xk = sorted[atomk].x;
      yk = sorted[atomk].y;
      zk = sorted[atomk].z;

      ijvdw = jvdw[i];
      imut = mut[i];
      kjvdw = jvdw[k];
      kmut = mut[k];

      for (int j = 0; j < WARP_SIZE; ++j) {
         bool incl = atomk > 0;
         real xr = xi - xk;
         real yr = yi - yk;
         real zr = zi - zk;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real r = REAL_SQRT(r2);
            real invr = REAL_RECIP(r);
            real rv = radmin[ijvdw * njvdw + kjvdw];
            real eps = epsilon[ijvdw * njvdw + kjvdw];
            real e, de;
            real vlambda = pair_vlambda(vlam, vcouple, imut, kmut);
            pair_lj_v3<do_g, true, 1>(r, invr, vlambda, 1, rv, eps, cut, off, e, de);
            if CONSTEXPR (do_e) {
               evtl += floatTo<ebuf_prec>(e);
               if CONSTEXPR (do_a) {
                  if (e != 0)
                     nevtl += 1;
               }
            }
            if CONSTEXPR (do_g) {
               real dedx, dedy, dedz;
               de = de * invr;
               dedx = de * xr;
               dedy = de * yr;
               dedz = de * zr;
               fix += dedx;
               fiy += dedy;
               fiz += dedz;
               fkx -= dedx;
               fky -= dedy;
               fkz -= dedz;
               if CONSTEXPR (do_v) {
                  vevtlxx += floatTo<vbuf_prec>(xr * dedx);
                  vevtlyx += floatTo<vbuf_prec>(yr * dedx);
                  vevtlzx += floatTo<vbuf_prec>(zr * dedx);
                  vevtlyy += floatTo<vbuf_prec>(yr * dedy);
                  vevtlzy += floatTo<vbuf_prec>(zr * dedy);
                  vevtlzz += floatTo<vbuf_prec>(zr * dedz);
               }
            }
         }

         xi = __shfl_sync(ALL_LANES, xi, ilane + 1);
         yi = __shfl_sync(ALL_LANES, yi, ilane + 1);
         zi = __shfl_sync(ALL_LANES, zi, ilane + 1);
         ijvdw = __shfl_sync(ALL_LANES, ijvdw, ilane + 1);
         imut = __shfl_sync(ALL_LANES, imut, ilane + 1);
         if CONSTEXPR (do_g) {
            fix = __shfl_sync(ALL_LANES, fix, ilane + 1);
            fiy = __shfl_sync(ALL_LANES, fiy, ilane + 1);
            fiz = __shfl_sync(ALL_LANES, fiz, ilane + 1);
         }
      }

      if CONSTEXPR (do_g) {
         atomic_add(fix, gx, i);
         atomic_add(fiy, gy, i);
         atomic_add(fiz, gz, i);
         atomic_add(fkx, gx, k);
         atomic_add(fky, gy, k);
         atomic_add(fkz, gz, k);
      }
   }

   if CONSTEXPR (do_a) {
      atomic_add(nevtl, nev, ithread);
   }
   if CONSTEXPR (do_e) {
      atomic_add(evtl, ev, ithread);
   }
   if CONSTEXPR (do_v) {
      atomic_add(vevtlxx, vevtlyx, vevtlzx, vevtlyy, vevtlzy, vevtlzz, vev, ithread);
   }
}

// special vdw14 interactions
template <class Ver>
__global__
void elj_cu2(CountBuffer restrict nebuf, EnergyBuffer restrict ebuf, VirialBuffer restrict vbuf,
   grad_prec* restrict gx, grad_prec* restrict gy, grad_prec* restrict gz, TINKER_IMAGE_PARAMS,
   real cut, real off, const real* restrict x, const real* restrict y, const real* restrict z, //
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

      int pos = jit * njvdw + jkt;
      real rv = radmin[pos];
      real eps = epsilon[pos];
      real rv4 = radmin4[pos];
      real eps4 = epsilon4[pos];

      real r2 = image2(xr, yr, zr);
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
   elj_cu1<Ver><<<ngrid, BLOCK_DIM, 0, g::s0>>>(st.n, TINKER_IMAGE_ARGS, nev, ev, vir_ev, devx,
      devy, devz, cut, off, st.si2.bit0, nvexclude, vexclude, vexclude_scale, st.x, st.y, st.z,
      st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst, njvdw, radmin, epsilon, jvdw, mut,
      vlam, vcouple);
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
