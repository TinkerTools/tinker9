#include "add.h"
#include "ff/energy.h"
#include "ff/hippo/edisp.h"
#include "ff/image.h"
#include "ff/pme.h"
#include "ff/spatial.h"
#include "ff/switch.h"
#include "launch.h"
#include "seq/bsplgen.h"
#include "seq/pair_disp.h"
#include "seq/triangle.h"
#include "tool/gpucard.h"

namespace tinker {
// ck.py Version 2.0.2
template <class Ver, class DTYP>
__global__
void edisp_cu1(int n, TINKER_IMAGE_PARAMS, CountBuffer restrict nd, EnergyBuffer restrict ed,
   VirialBuffer restrict vd, grad_prec* restrict gx, grad_prec* restrict gy, grad_prec* restrict gz,
   real cut, real off, const unsigned* restrict dinfo, int nexclude,
   const int (*restrict exclude)[2], const real* restrict exclude_scale, const real* restrict x,
   const real* restrict y, const real* restrict z, const Spatial::SortedAtom* restrict sorted,
   int nakpl, const int* restrict iakpl, int niak, const int* restrict iak, const int* restrict lst,
   const real* restrict csix, const real* restrict adisp, real aewald)
{
   constexpr bool do_a = Ver::a;
   constexpr bool do_e = Ver::e;
   constexpr bool do_v = Ver::v;
   constexpr bool do_g = Ver::g;
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);

   int ndtl;
   if CONSTEXPR (do_a) {
      ndtl = 0;
   }
   using ebuf_prec = EnergyBufferTraits::type;
   ebuf_prec edtl;
   if CONSTEXPR (do_e) {
      edtl = 0;
   }
   using vbuf_prec = VirialBufferTraits::type;
   vbuf_prec vdtlxx, vdtlyx, vdtlzx, vdtlyy, vdtlzy, vdtlzz;
   if CONSTEXPR (do_v) {
      vdtlxx = 0;
      vdtlyx = 0;
      vdtlzx = 0;
      vdtlyy = 0;
      vdtlzy = 0;
      vdtlzz = 0;
   }
   real xi;
   real yi;
   real zi;
   real xk;
   real yk;
   real zk;
   real gxi;
   real gyi;
   real gzi;
   real gxk;
   real gyk;
   real gzk;
   real ci;
   real ai;
   real ck;
   real ak;

   //* /
   for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
      if CONSTEXPR (do_g) {
         gxi = 0;
         gyi = 0;
         gzi = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
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
      ci = csix[i];
      ai = adisp[i];
      ck = csix[k];
      ak = adisp[k];

      constexpr bool incl = true;
      real xr = xi - xk;
      real yr = yi - yk;
      real zr = zi - zk;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         real r = REAL_SQRT(r2);
         real rr1 = REAL_RECIP(r);
         real e, de;
         pair_disp<do_g, DTYP, 0>(r, r2, rr1, scalea, aewald, ci, ai, ck, ak, cut, off, e, de);
         if CONSTEXPR (do_e) {
            edtl += floatTo<ebuf_prec>(e);
            if CONSTEXPR (do_a) {
               if (scalea != 0 and e != 0)
                  ndtl += 1;
            }
         }
         if CONSTEXPR (do_g) {
            real dedx, dedy, dedz;
            de *= rr1;
            dedx = de * xr;
            dedy = de * yr;
            dedz = de * zr;
            gxi += dedx;
            gyi += dedy;
            gzi += dedz;
            gxk -= dedx;
            gyk -= dedy;
            gzk -= dedz;
            if CONSTEXPR (do_v) {
               vdtlxx += floatTo<vbuf_prec>(xr * dedx);
               vdtlyx += floatTo<vbuf_prec>(yr * dedx);
               vdtlzx += floatTo<vbuf_prec>(zr * dedx);
               vdtlyy += floatTo<vbuf_prec>(yr * dedy);
               vdtlzy += floatTo<vbuf_prec>(zr * dedy);
               vdtlzz += floatTo<vbuf_prec>(zr * dedz);
            }
         }
      } // end if (include)

      if CONSTEXPR (do_g) {
         atomic_add(gxi, gx, i);
         atomic_add(gyi, gy, i);
         atomic_add(gzi, gz, i);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
      }
   }
   // */

   for (int iw = iwarp; iw < nakpl; iw += nwarp) {
      if CONSTEXPR (do_g) {
         gxi = 0;
         gyi = 0;
         gzi = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
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

      ci = csix[i];
      ai = adisp[i];
      ck = csix[k];
      ak = adisp[k];

      unsigned int dinfo0 = dinfo[iw * WARP_SIZE + ilane];
      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         bool incl = iid < kid and kid < n;
         int srcmask = 1 << srclane;
         incl = incl and (dinfo0 & srcmask) == 0;
         real xr = xi - xk;
         real yr = yi - yk;
         real zr = zi - zk;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real r = REAL_SQRT(r2);
            real rr1 = REAL_RECIP(r);
            real e, de;
            pair_disp<do_g, DTYP, 1>(r, r2, rr1, 1, aewald, ci, ai, ck, ak, cut, off, e, de);
            if CONSTEXPR (do_e) {
               edtl += floatTo<ebuf_prec>(e);
               if CONSTEXPR (do_a) {
                  if (e != 0)
                     ndtl += 1;
               }
            }
            if CONSTEXPR (do_g) {
               real dedx, dedy, dedz;
               de *= rr1;
               dedx = de * xr;
               dedy = de * yr;
               dedz = de * zr;
               gxi += dedx;
               gyi += dedy;
               gzi += dedz;
               gxk -= dedx;
               gyk -= dedy;
               gzk -= dedz;
               if CONSTEXPR (do_v) {
                  vdtlxx += floatTo<vbuf_prec>(xr * dedx);
                  vdtlyx += floatTo<vbuf_prec>(yr * dedx);
                  vdtlzx += floatTo<vbuf_prec>(zr * dedx);
                  vdtlyy += floatTo<vbuf_prec>(yr * dedy);
                  vdtlzy += floatTo<vbuf_prec>(zr * dedy);
                  vdtlzz += floatTo<vbuf_prec>(zr * dedz);
               }
            }
         } // end if (include)

         iid = __shfl_sync(ALL_LANES, iid, ilane + 1);
         xi = __shfl_sync(ALL_LANES, xi, ilane + 1);
         yi = __shfl_sync(ALL_LANES, yi, ilane + 1);
         zi = __shfl_sync(ALL_LANES, zi, ilane + 1);
         ci = __shfl_sync(ALL_LANES, ci, ilane + 1);
         ai = __shfl_sync(ALL_LANES, ai, ilane + 1);
         if CONSTEXPR (do_g) {
            gxi = __shfl_sync(ALL_LANES, gxi, ilane + 1);
            gyi = __shfl_sync(ALL_LANES, gyi, ilane + 1);
            gzi = __shfl_sync(ALL_LANES, gzi, ilane + 1);
         }
      }

      if CONSTEXPR (do_g) {
         atomic_add(gxi, gx, i);
         atomic_add(gyi, gy, i);
         atomic_add(gzi, gz, i);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
      }
   }

   for (int iw = iwarp; iw < niak; iw += nwarp) {
      if CONSTEXPR (do_g) {
         gxi = 0;
         gyi = 0;
         gzi = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
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

      ci = csix[i];
      ai = adisp[i];
      ck = csix[k];
      ak = adisp[k];

      for (int j = 0; j < WARP_SIZE; ++j) {
         bool incl = atomk > 0;
         real xr = xi - xk;
         real yr = yi - yk;
         real zr = zi - zk;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real r = REAL_SQRT(r2);
            real rr1 = REAL_RECIP(r);
            real e, de;
            pair_disp<do_g, DTYP, 1>(r, r2, rr1, 1, aewald, ci, ai, ck, ak, cut, off, e, de);
            if CONSTEXPR (do_e) {
               edtl += floatTo<ebuf_prec>(e);
               if CONSTEXPR (do_a) {
                  if (e != 0)
                     ndtl += 1;
               }
            }
            if CONSTEXPR (do_g) {
               real dedx, dedy, dedz;
               de *= rr1;
               dedx = de * xr;
               dedy = de * yr;
               dedz = de * zr;
               gxi += dedx;
               gyi += dedy;
               gzi += dedz;
               gxk -= dedx;
               gyk -= dedy;
               gzk -= dedz;
               if CONSTEXPR (do_v) {
                  vdtlxx += floatTo<vbuf_prec>(xr * dedx);
                  vdtlyx += floatTo<vbuf_prec>(yr * dedx);
                  vdtlzx += floatTo<vbuf_prec>(zr * dedx);
                  vdtlyy += floatTo<vbuf_prec>(yr * dedy);
                  vdtlzy += floatTo<vbuf_prec>(zr * dedy);
                  vdtlzz += floatTo<vbuf_prec>(zr * dedz);
               }
            }
         } // end if (include)

         xi = __shfl_sync(ALL_LANES, xi, ilane + 1);
         yi = __shfl_sync(ALL_LANES, yi, ilane + 1);
         zi = __shfl_sync(ALL_LANES, zi, ilane + 1);
         ci = __shfl_sync(ALL_LANES, ci, ilane + 1);
         ai = __shfl_sync(ALL_LANES, ai, ilane + 1);
         if CONSTEXPR (do_g) {
            gxi = __shfl_sync(ALL_LANES, gxi, ilane + 1);
            gyi = __shfl_sync(ALL_LANES, gyi, ilane + 1);
            gzi = __shfl_sync(ALL_LANES, gzi, ilane + 1);
         }
      }

      if CONSTEXPR (do_g) {
         atomic_add(gxi, gx, i);
         atomic_add(gyi, gy, i);
         atomic_add(gzi, gz, i);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
      }
   }

   if CONSTEXPR (do_a) {
      atomic_add(ndtl, nd, ithread);
   }
   if CONSTEXPR (do_e) {
      atomic_add(edtl, ed, ithread);
   }
   if CONSTEXPR (do_v) {
      atomic_add(vdtlxx, vdtlyx, vdtlzx, vdtlyy, vdtlzy, vdtlzz, vd, ithread);
   }
}

template <class Ver, class DTYP>
void edisp_cu()
{
   const auto& st = *dspspatial_v2_unit;
   real cut, off;
   if CONSTEXPR (eq<DTYP, DEWALD>()) {
      off = switchOff(Switch::DEWALD);
      cut = off; // not used
   } else {
      off = switchOff(Switch::DISP);
      cut = switchCut(Switch::DISP);
   }

   real aewald = 0;
   if CONSTEXPR (eq<DTYP, DEWALD>()) {
      PMEUnit pu = dpme_unit;
      aewald = pu->aewald;
   }

   int ngrid = gpuGridSize(BLOCK_DIM);
   edisp_cu1<Ver, DTYP><<<ngrid, BLOCK_DIM, 0, g::s0>>>(st.n, TINKER_IMAGE_ARGS, ndisp, edsp,
      vir_edsp, dedspx, dedspy, dedspz, cut, off, st.si1.bit0, ndspexclude, dspexclude,
      dspexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst,
      csix, adisp, aewald);
}

void edisp_ewald_real_cu(int vers)
{
   if (vers == calc::v0)
      edisp_cu<calc::V0, DEWALD>();
   else if (vers == calc::v1)
      edisp_cu<calc::V1, DEWALD>();
   else if (vers == calc::v3)
      edisp_cu<calc::V3, DEWALD>();
   else if (vers == calc::v4)
      edisp_cu<calc::V4, DEWALD>();
   else if (vers == calc::v5)
      edisp_cu<calc::V5, DEWALD>();
   else if (vers == calc::v6)
      edisp_cu<calc::V6, DEWALD>();
}

template <class Ver, int bsorder>
__global__
void edisp_cu3(CountBuffer restrict ndisp, EnergyBuffer restrict edsp, const real* restrict csix,
   real aewald, int n, int nfft1, int nfft2, int nfft3, const real* restrict x,
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
      real icsix = csix[ii];
      if (icsix == 0)
         continue;

      // self energy
      if CONSTEXPR (do_e) {
         real fs = aewald * aewald;
         fs *= fs * fs;
         fs /= 12;
         real e = fs * icsix * icsix;
         atomic_add(e, edsp, ithread);
         if CONSTEXPR (do_a) {
            atomic_add(1, ndisp, ithread);
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

         real fi = csix[ii];
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
                  real t1 = thetai1[4 * ix];
                  real dt1 = nfft1 * thetai1[1 + 4 * ix];
                  real term = qgrid[2 * index];
                  de1 += 2 * term * dt1 * t2 * t3;
                  de2 += 2 * term * dt2 * t1 * t3;
                  de3 += 2 * term * dt3 * t1 * t2;
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
void edisp_cu4()
{
   const auto& st = *dpme_unit;
   real aewald = st.aewald;
   int nfft1 = st.nfft1;
   int nfft2 = st.nfft2;
   int nfft3 = st.nfft3;

   assert(st.bsorder == 4);
   auto ker = edisp_cu3<Ver, 4>;
   launch_k2b(g::s0, PME_BLOCKDIM, n, ker, //
      ndisp, edsp, csix, aewald, n, nfft1, nfft2, nfft3, x, y, z, st.qgrid, recipa, recipb, recipc,
      dedspx, dedspy, dedspz);
}

void edisp_ewald_recip_self_cu(int vers)
{
   if (vers == calc::v0)
      edisp_cu4<calc::V0>();
   else if (vers == calc::v1)
      edisp_cu4<calc::V1>();
   else if (vers == calc::v3)
      edisp_cu4<calc::V3>();
   else if (vers == calc::v4)
      edisp_cu4<calc::V4>();
   else if (vers == calc::v5)
      edisp_cu4<calc::V5>();
   else if (vers == calc::v6)
      edisp_cu4<calc::V6>();
}

void edisp_nonewald_cu(int vers)
{
   if (vers == calc::v0)
      edisp_cu<calc::V0, NON_EWALD_TAPER>();
   else if (vers == calc::v1)
      edisp_cu<calc::V1, NON_EWALD_TAPER>();
   else if (vers == calc::v3)
      edisp_cu<calc::V3, NON_EWALD_TAPER>();
   else if (vers == calc::v4)
      edisp_cu<calc::V4, NON_EWALD_TAPER>();
   else if (vers == calc::v5)
      edisp_cu<calc::V5, NON_EWALD_TAPER>();
   else if (vers == calc::v6)
      edisp_cu<calc::V6, NON_EWALD_TAPER>();
}
}
