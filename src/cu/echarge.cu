#include "ff/energy.h"
#include "ff/image.h"
#include "ff/pchg/echarge.h"
#include "ff/pme.h"
#include "ff/spatial.h"
#include "ff/switch.h"
#include "launch.h"
#include "seq/bsplgen.h"
#include "seq/pair_charge.h"
#include "triangle.h"

namespace tinker {
// ck.py Version 2.0.2
template <class Ver, class ETYP>
__global__
void echarge_cu1(int n, TINKER_IMAGE_PARAMS, CountBuffer restrict nec, EnergyBuffer restrict ec,
   VirialBuffer restrict vec, grad_prec* restrict gx, grad_prec* restrict gy,
   grad_prec* restrict gz, real cut, real off, const unsigned* restrict info, int nexclude,
   const int (*restrict exclude)[2], const real* restrict exclude_scale, const real* restrict x,
   const real* restrict y, const real* restrict z, const Spatial::SortedAtom* restrict sorted,
   int nakpl, const int* restrict iakpl, int niak, const int* restrict iak, const int* restrict lst,
   real ebuffer, real f, real aewald, const real* restrict chg)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);

   int nectl;
   if CONSTEXPR (do_a) {
      nectl = 0;
   }
   using ebuf_prec = EnergyBufferTraits::type;
   ebuf_prec ectl;
   if CONSTEXPR (do_e) {
      ectl = 0;
   }
   using vbuf_prec = VirialBufferTraits::type;
   vbuf_prec vectlxx, vectlyx, vectlzx, vectlyy, vectlzy, vectlzz;
   if CONSTEXPR (do_v) {
      vectlxx = 0;
      vectlyx = 0;
      vectlzx = 0;
      vectlyy = 0;
      vectlzy = 0;
      vectlzz = 0;
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
   real ichg;
   real kchg;

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
      ichg = chg[i];
      kchg = chg[k];

      constexpr bool incl = true;
      real xr = xi - xk;
      real yr = yi - yk;
      real zr = zi - zk;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         real r = REAL_SQRT(r2);
         real invr = REAL_RECIP(r);
         real e, de;
         pair_chg_v3<do_g, ETYP, 0>(r, scalea, ichg, kchg, ebuffer, f, aewald, cut, off, e, de);
         if CONSTEXPR (do_e) {
            ectl += floatTo<ebuf_prec>(e);
            if CONSTEXPR (do_a) {
               if (scalea != 0 and e != 0)
                  nectl += 1;
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
               vectlxx += floatTo<vbuf_prec>(xr * dedx);
               vectlyx += floatTo<vbuf_prec>(yr * dedx);
               vectlzx += floatTo<vbuf_prec>(zr * dedx);
               vectlyy += floatTo<vbuf_prec>(yr * dedy);
               vectlzy += floatTo<vbuf_prec>(zr * dedy);
               vectlzz += floatTo<vbuf_prec>(zr * dedz);
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

      ichg = chg[i];
      kchg = chg[k];

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
            real e, de;
            pair_chg_v3<do_g, ETYP, 1>(r, 1, ichg, kchg, ebuffer, f, aewald, cut, off, e, de);
            if CONSTEXPR (do_e) {
               ectl += floatTo<ebuf_prec>(e);
               if CONSTEXPR (do_a) {
                  if (e != 0)
                     nectl += 1;
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
                  vectlxx += floatTo<vbuf_prec>(xr * dedx);
                  vectlyx += floatTo<vbuf_prec>(yr * dedx);
                  vectlzx += floatTo<vbuf_prec>(zr * dedx);
                  vectlyy += floatTo<vbuf_prec>(yr * dedy);
                  vectlzy += floatTo<vbuf_prec>(zr * dedy);
                  vectlzz += floatTo<vbuf_prec>(zr * dedz);
               }
            }
         }

         iid = __shfl_sync(ALL_LANES, iid, ilane + 1);
         xi = __shfl_sync(ALL_LANES, xi, ilane + 1);
         yi = __shfl_sync(ALL_LANES, yi, ilane + 1);
         zi = __shfl_sync(ALL_LANES, zi, ilane + 1);
         ichg = __shfl_sync(ALL_LANES, ichg, ilane + 1);
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

      ichg = chg[i];
      kchg = chg[k];

      for (int j = 0; j < WARP_SIZE; ++j) {
         bool incl = atomk > 0;
         real xr = xi - xk;
         real yr = yi - yk;
         real zr = zi - zk;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real r = REAL_SQRT(r2);
            real invr = REAL_RECIP(r);
            real e, de;
            pair_chg_v3<do_g, ETYP, 1>(r, 1, ichg, kchg, ebuffer, f, aewald, cut, off, e, de);
            if CONSTEXPR (do_e) {
               ectl += floatTo<ebuf_prec>(e);
               if CONSTEXPR (do_a) {
                  if (e != 0)
                     nectl += 1;
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
                  vectlxx += floatTo<vbuf_prec>(xr * dedx);
                  vectlyx += floatTo<vbuf_prec>(yr * dedx);
                  vectlzx += floatTo<vbuf_prec>(zr * dedx);
                  vectlyy += floatTo<vbuf_prec>(yr * dedy);
                  vectlzy += floatTo<vbuf_prec>(zr * dedy);
                  vectlzz += floatTo<vbuf_prec>(zr * dedz);
               }
            }
         }

         xi = __shfl_sync(ALL_LANES, xi, ilane + 1);
         yi = __shfl_sync(ALL_LANES, yi, ilane + 1);
         zi = __shfl_sync(ALL_LANES, zi, ilane + 1);
         ichg = __shfl_sync(ALL_LANES, ichg, ilane + 1);
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
      atomic_add(nectl, nec, ithread);
   }
   if CONSTEXPR (do_e) {
      atomic_add(ectl, ec, ithread);
   }
   if CONSTEXPR (do_v) {
      atomic_add(vectlxx, vectlyx, vectlzx, vectlyy, vectlzy, vectlzz, vec, ithread);
   }
}

template <class Ver, class ETYP>
static void echarge_cu()
{
   const auto& st = *cspatial_v2_unit;
   real cut, off;
   if CONSTEXPR (eq<ETYP, EWALD>()) {
      off = switchOff(Switch::EWALD);
      cut = off; // not used
   } else {
      off = switchOff(Switch::CHARGE);
      cut = switchCut(Switch::CHARGE);
   }

   const real f = electric / dielec;
   real aewald = 0;
   if CONSTEXPR (eq<ETYP, EWALD>()) {
      PMEUnit pu = epme_unit;
      aewald = pu->aewald;
   }

   int ngrid = gpuGridSize(BLOCK_DIM);
   auto ker1 = echarge_cu1<Ver, ETYP>;
   ker1<<<ngrid, BLOCK_DIM, 0, g::s0>>>(st.n, TINKER_IMAGE_ARGS, nec, ec, vir_ec, decx, decy, decz,
      cut, off, st.si1.bit0, ncexclude, cexclude, cexclude_scale, st.x, st.y, st.z, st.sorted,
      st.nakpl, st.iakpl, st.niak, st.iak, st.lst, ebuffer, f, aewald, pchg);
}

void echargeNonEwald_cu(int vers)
{
   if (vers == calc::v0)
      echarge_cu<calc::V0, NON_EWALD_TAPER>();
   else if (vers == calc::v1)
      echarge_cu<calc::V1, NON_EWALD_TAPER>();
   else if (vers == calc::v3)
      echarge_cu<calc::V3, NON_EWALD_TAPER>();
   else if (vers == calc::v4)
      echarge_cu<calc::V4, NON_EWALD_TAPER>();
   else if (vers == calc::v5)
      echarge_cu<calc::V5, NON_EWALD_TAPER>();
   else if (vers == calc::v6)
      echarge_cu<calc::V6, NON_EWALD_TAPER>();
}

void echargeEwaldReal_cu(int vers)
{
   if (vers == calc::v0)
      echarge_cu<calc::V0, EWALD>();
   else if (vers == calc::v1)
      echarge_cu<calc::V1, EWALD>();
   else if (vers == calc::v3)
      echarge_cu<calc::V3, EWALD>();
   else if (vers == calc::v4)
      echarge_cu<calc::V4, EWALD>();
   else if (vers == calc::v5)
      echarge_cu<calc::V5, EWALD>();
   else if (vers == calc::v6)
      echarge_cu<calc::V6, EWALD>();
}

template <class Ver, int bsorder>
__global__
void echarge_cu3(CountBuffer restrict nec, EnergyBuffer restrict ec, const real* restrict pchg,
   real f, real aewald, int n, int nfft1, int nfft2, int nfft3, const real* restrict x,
   const real* restrict y, const real* restrict z, const real* restrict qgrid, real3 reca,
   real3 recb, real3 recc, grad_prec* restrict gx, grad_prec* restrict gy, grad_prec* restrict gz)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;

   real thetai1[4 * 5];
   real thetai2[4 * 5];
   real thetai3[4 * 5];
   __shared__ real sharedarray[5 * 5 * PME_BLOCKDIM];
   real* restrict array = &sharedarray[5 * 5 * threadIdx.x];

   for (int ii = ithread; ii < n; ii += blockDim.x * gridDim.x) {
      real chgi = pchg[ii];
      if (chgi == 0)
         continue;

      // self energy, tinfoil
      if CONSTEXPR (do_e) {
         real fs = -f * aewald * REAL_RECIP(sqrtpi);
         real e = fs * chgi * chgi;
         atomic_add(e, ec, ithread);
         if CONSTEXPR (do_a) {
            atomic_add(1, nec, ithread);
         }
      }

      // recip gradient
      if CONSTEXPR (do_g) {
         real xi = x[ii];
         real yi = y[ii];
         real zi = z[ii];

         real w1 = xi * reca.x + yi * reca.y + zi * reca.z;
         w1 = w1 + 0.5f - REAL_FLOOR(w1 + 0.5f);
         real fr1 = nfft1 * w1;
         int igrid1 = REAL_FLOOR(fr1);
         w1 = fr1 - igrid1;

         real w2 = xi * recb.x + yi * recb.y + zi * recb.z;
         w2 = w2 + 0.5f - REAL_FLOOR(w2 + 0.5f);
         real fr2 = nfft2 * w2;
         int igrid2 = REAL_FLOOR(fr2);
         w2 = fr2 - igrid2;

         real w3 = xi * recc.x + yi * recc.y + zi * recc.z;
         w3 = w3 + 0.5f - REAL_FLOOR(w3 + 0.5f);
         real fr3 = nfft3 * w3;
         int igrid3 = REAL_FLOOR(fr3);
         w3 = fr3 - igrid3;

         igrid1 = igrid1 - bsorder + 1;
         igrid2 = igrid2 - bsorder + 1;
         igrid3 = igrid3 - bsorder + 1;
         igrid1 += (igrid1 < 0 ? nfft1 : 0);
         igrid2 += (igrid2 < 0 ? nfft2 : 0);
         igrid3 += (igrid3 < 0 ? nfft3 : 0);

         bsplgen<2, bsorder>(w1, thetai1, array);
         bsplgen<2, bsorder>(w2, thetai2, array);
         bsplgen<2, bsorder>(w3, thetai3, array);

         real fi = f * chgi;
         real de1 = 0, de2 = 0, de3 = 0;
         for (int iz = 0; iz < bsorder; ++iz) {
            int zbase = igrid3 + iz;
            zbase -= (zbase >= nfft3 ? nfft3 : 0);
            zbase *= (nfft1 * nfft2);
            real t3 = thetai3[4 * iz];
            real dt3 = nfft3 * thetai3[1 + 4 * iz];
            for (int iy = 0; iy < bsorder; ++iy) {
               int ybase = igrid2 + iy;
               ybase -= (ybase >= nfft2 ? nfft2 : 0);
               ybase *= nfft1;
               real t2 = thetai2[4 * iy];
               real dt2 = nfft2 * thetai2[1 + 4 * iy];
               for (int ix = 0; ix < bsorder; ++ix) {
                  int xbase = igrid1 + ix;
                  xbase -= (xbase >= nfft1 ? nfft1 : 0);
                  int index = xbase + ybase + zbase;
                  real term = qgrid[2 * index];
                  real t1 = thetai1[4 * ix];
                  real dt1 = nfft1 * thetai1[1 + 4 * ix];
                  de1 += term * dt1 * t2 * t3;
                  de2 += term * dt2 * t1 * t3;
                  de3 += term * dt3 * t1 * t2;
               }
            }
         } // end for (iz)

         real frcx = fi * (reca.x * de1 + recb.x * de2 + recc.x * de3);
         real frcy = fi * (reca.y * de1 + recb.y * de2 + recc.y * de3);
         real frcz = fi * (reca.z * de1 + recb.z * de2 + recc.z * de3);
         atomic_add(frcx, gx, ii);
         atomic_add(frcy, gy, ii);
         atomic_add(frcz, gz, ii);
      }
   }
}

template <class Ver>
static void echargeFphiSelf_cu()
{
   real f = electric / dielec;
   const auto& st = *epme_unit;
   real aewald = st.aewald;
   int nfft1 = st.nfft1;
   int nfft2 = st.nfft2;
   int nfft3 = st.nfft3;

   auto stream = g::s0;
   if (use_pme_stream)
      stream = g::spme;
   if (st.bsorder == 5) {
      auto ker = echarge_cu3<Ver, 5>;
      launch_k2b(stream, PME_BLOCKDIM, n, ker, //
         nec, ec, pchg, f, aewald, n,          //
         nfft1, nfft2, nfft3, x, y, z, st.qgrid, recipa, recipb, recipc, decx, decy, decz);
   } else if (st.bsorder == 4) {
      auto ker = echarge_cu3<Ver, 4>;
      launch_k2b(stream, PME_BLOCKDIM, n, ker, //
         nec, ec, pchg, f, aewald, n,          //
         nfft1, nfft2, nfft3, x, y, z, st.qgrid, recipa, recipb, recipc, decx, decy, decz);
   }
}

void echargeEwaldFphiSelf_cu(int vers)
{
   if (vers == calc::v0)
      echargeFphiSelf_cu<calc::V0>();
   else if (vers == calc::v1)
      echargeFphiSelf_cu<calc::V1>();
   else if (vers == calc::v3)
      echargeFphiSelf_cu<calc::V3>();
   else if (vers == calc::v4)
      echargeFphiSelf_cu<calc::V4>();
   else if (vers == calc::v5)
      echargeFphiSelf_cu<calc::V5>();
   else if (vers == calc::v6)
      echargeFphiSelf_cu<calc::V6>();
}
}
