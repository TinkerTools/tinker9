#include "add.h"
#include "edisp.h"
#include "glob.spatial.h"
#include "image.h"
#include "launch.h"
#include "md.h"
#include "pmestuf.h"
#include "seq_bsplgen.h"
#include "seq_pair_disp.h"
#include "seq_switch.h"
#include "switch.h"


namespace tinker {
#define EDISP_ARGS                                                             \
   size_t bufsize, count_buffer restrict ndisp, energy_buffer restrict edsp,   \
      virial_buffer restrict vir_edsp, grad_prec *restrict gx,                 \
      grad_prec *restrict gy, grad_prec *restrict gz, TINKER_IMAGE_PARAMS,     \
      real cut, real off, const real *restrict csix,                           \
      const real *restrict adisp


template <class Ver, class DTYP>
__global__
void edisp_cu1(EDISP_ARGS, const Spatial::SortedAtom* restrict sorted, int niak,
               const int* restrict iak, const int* restrict lst, int n,
               real aewald)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);
   const int offset = ithread & (bufsize - 1);


   MAYBE_UNUSED int ctl;
   MAYBE_UNUSED e_prec etl;
   MAYBE_UNUSED v_prec vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz;
   MAYBE_UNUSED real gxi, gyi, gzi, gxk, gyk, gzk;


   const real cut2 = cut * cut;
   const real off2 = off * off;
   for (int iw = iwarp; iw < niak; iw += nwarp) {
      if CONSTEXPR (do_a)
         ctl = 0;
      if CONSTEXPR (do_e)
         etl = 0;
      if CONSTEXPR (do_g) {
         gxi = 0;
         gyi = 0;
         gzi = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
      }
      if CONSTEXPR (do_v) {
         vtlxx = 0;
         vtlyx = 0;
         vtlzx = 0;
         vtlyy = 0;
         vtlzy = 0;
         vtlzz = 0;
      }


      int atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
      real xi = sorted[atomi].x;
      real yi = sorted[atomi].y;
      real zi = sorted[atomi].z;
      int i = sorted[atomi].unsorted;
      real ci = csix[i];
      real ai = adisp[i];


      int shatomk = lst[iw * WARP_SIZE + ilane];
      real shx = sorted[shatomk].x;
      real shy = sorted[shatomk].y;
      real shz = sorted[shatomk].z;
      int shk = sorted[shatomk].unsorted;
      real shck = csix[shk];
      real shak = adisp[shk];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
         real xr = xi - __shfl_sync(ALL_LANES, shx, srclane);
         real yr = yi - __shfl_sync(ALL_LANES, shy, srclane);
         real zr = zi - __shfl_sync(ALL_LANES, shz, srclane);
         real ck = __shfl_sync(ALL_LANES, shck, srclane);
         real ak = __shfl_sync(ALL_LANES, shak, srclane);


         MAYBE_UNUSED real dedx = 0, dedy = 0, dedz = 0;
         real r2 = image2(xr, yr, zr);
         if (atomi < atomk && r2 <= off2) {
            real r = REAL_SQRT(r2);


            MAYBE_UNUSED e_prec e, de;
            if CONSTEXPR (eq<DTYP, DEWALD>()) {
               pair_disp<do_g, DEWALD>(r, r2, 1, aewald, ci, ai, ck, ak, e, de);
            } else if CONSTEXPR (eq<DTYP, NON_EWALD_TAPER>()) {
               pair_disp<do_g, NON_EWALD_TAPER>(r, r2, 1, 0, ci, ai, ck, ak, e,
                                                de);
               if (r2 > cut2) {
                  real taper, dtaper;
                  switch_taper5<do_g>(r, cut, off, taper, dtaper);
                  if CONSTEXPR (do_g)
                     de = e * dtaper + de * taper;
                  if CONSTEXPR (do_e)
                     e *= taper;
               }
            }


            if CONSTEXPR (do_a)
               if (e != 0)
                  ctl += 1;
            if CONSTEXPR (do_e)
               etl += e;
            if CONSTEXPR (do_g) {
               de *= REAL_RECIP(r);
               dedx = de * xr;
               dedy = de * yr;
               dedz = de * zr;
               if CONSTEXPR (do_v) {
                  vtlxx += xr * dedx;
                  vtlyx += yr * dedx;
                  vtlzx += zr * dedx;
                  vtlyy += yr * dedy;
                  vtlzy += zr * dedy;
                  vtlzz += zr * dedz;
               }
            }
         } // end if (include)


         if CONSTEXPR (do_g) {
            int dstlane = (ilane + WARP_SIZE - j) & (WARP_SIZE - 1);
            gxi += dedx;
            gyi += dedy;
            gzi += dedz;
            gxk -= __shfl_sync(ALL_LANES, dedx, dstlane);
            gyk -= __shfl_sync(ALL_LANES, dedy, dstlane);
            gzk -= __shfl_sync(ALL_LANES, dedz, dstlane);
         }
      }


      if CONSTEXPR (do_a)
         atomic_add(ctl, ndisp, offset);
      if CONSTEXPR (do_e)
         atomic_add(etl, edsp, offset);
      if CONSTEXPR (do_g) {
         atomic_add(gxi, gx, i);
         atomic_add(gyi, gy, i);
         atomic_add(gzi, gz, i);
         atomic_add(gxk, gx, shk);
         atomic_add(gyk, gy, shk);
         atomic_add(gzk, gz, shk);
      }
      if CONSTEXPR (do_v)
         atomic_add(vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz, vir_edsp, offset);
   } // end for (iw)
}


template <class Ver, class DTYP>
__global__
void edisp_cu2(EDISP_ARGS, const real* restrict x, const real* restrict y,
               const real* restrict z, int ndspexclude,
               const int (*restrict dspexclude)[2],
               const real* restrict dspexclude_scale)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   const real cut2 = cut * cut;
   const real off2 = off * off;
   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < ndspexclude;
        ii += blockDim.x * gridDim.x) {
      int offset = ii & (bufsize - 1);


      int i = dspexclude[ii][0];
      int k = dspexclude[ii][1];
      real dspscale = dspexclude_scale[ii];


      real ci = csix[i];
      real ai = adisp[i];
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];


      real ck = csix[k];
      real ak = adisp[k];
      real xr = xi - x[k];
      real yr = yi - y[k];
      real zr = zi - z[k];


      real r2 = image2(xr, yr, zr);
      if (r2 <= off2) {
         real r = REAL_SQRT(r2);


         MAYBE_UNUSED e_prec e, de;
         pair_disp<do_g, NON_EWALD>(r, r2, dspscale, 0, ci, ai, ck, ak, e, de);
         if CONSTEXPR (eq<DTYP, NON_EWALD_TAPER>()) {
            if (r2 > cut2) {
               real taper, dtaper;
               switch_taper5<do_g>(r, cut, off, taper, dtaper);
               if CONSTEXPR (do_g)
                  de = e * dtaper + de * taper;
               if CONSTEXPR (do_e)
                  e = e * taper;
            }
         }


         if CONSTEXPR (do_a)
            if (dspscale == -1 && e != 0)
               atomic_add(-1, ndisp, offset);
         if CONSTEXPR (do_e)
            atomic_add(e, edsp, offset);
         if CONSTEXPR (do_g) {
            de *= REAL_RECIP(r);
            real dedx = de * xr;
            real dedy = de * yr;
            real dedz = de * zr;
            atomic_add(dedx, gx, i);
            atomic_add(dedy, gy, i);
            atomic_add(dedz, gz, i);
            atomic_add(-dedx, gx, k);
            atomic_add(-dedy, gy, k);
            atomic_add(-dedz, gz, k);
            if CONSTEXPR (do_v) {
               real vxx = xr * dedx;
               real vyx = yr * dedx;
               real vzx = zr * dedx;
               real vyy = yr * dedy;
               real vzy = zr * dedy;
               real vzz = zr * dedz;
               atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_edsp, offset);
            }
         }
      }
   }
}


template <class Ver, class DTYP>
void edisp_cu()
{
   const auto& st = *dspspatial_unit;
   real cut, off;
   if CONSTEXPR (eq<DTYP, DEWALD>()) {
      off = switch_off(switch_dewald);
      // cut = off; // not used
   } else {
      off = switch_off(switch_disp);
      cut = switch_cut(switch_disp);
   }
   size_t bufsize = buffer_size();


   real aewald = 0;
   if CONSTEXPR (eq<DTYP, DEWALD>()) {
      PMEUnit pu = dpme_unit;
      aewald = pu->aewald;
   }


   if CONSTEXPR (eq<DTYP, DEWALD>()) {
      if (st.niak > 0) {
         auto ker1 = edisp_cu1<Ver, DEWALD>;
         launch_k1s(nonblk, WARP_SIZE * st.niak, ker1, //
                    bufsize, ndisp, edsp, vir_edsp, dedspx, dedspy, dedspz,
                    TINKER_IMAGE_ARGS, 0, off, csix, adisp, //
                    st.sorted, st.niak, st.iak, st.lst, n, aewald);
      }
      if (ndspexclude > 0) {
         auto ker2 = edisp_cu2<Ver, NON_EWALD>;
         launch_k1s(nonblk, ndspexclude, ker2, //
                    bufsize, ndisp, edsp, vir_edsp, dedspx, dedspy, dedspz,
                    TINKER_IMAGE_ARGS, 0, off, csix, adisp, //
                    x, y, z, ndspexclude, dspexclude, dspexclude_scale);
      }
   } else if CONSTEXPR (eq<DTYP, NON_EWALD_TAPER>()) {
      if (st.niak > 0) {
         auto ker1 = edisp_cu1<Ver, NON_EWALD_TAPER>;
         launch_k1s(nonblk, WARP_SIZE * st.niak, ker1, //
                    bufsize, ndisp, edsp, vir_edsp, dedspx, dedspy, dedspz,
                    TINKER_IMAGE_ARGS, cut, off, csix, adisp, //
                    st.sorted, st.niak, st.iak, st.lst, n, aewald);
      }
      if (ndspexclude > 0) {
         auto ker2 = edisp_cu2<Ver, NON_EWALD_TAPER>;
         launch_k1s(nonblk, ndspexclude, ker2, //
                    bufsize, ndisp, edsp, vir_edsp, dedspx, dedspy, dedspz,
                    TINKER_IMAGE_ARGS, cut, off, csix, adisp, //
                    x, y, z, ndspexclude, dspexclude, dspexclude_scale);
      }
   }
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
void edisp_cu3(size_t bufsize, count_buffer restrict ndisp,
               energy_buffer restrict edsp, const real* restrict csix,
               real aewald, int n, int nfft1, int nfft2, int nfft3,
               const real* restrict x, const real* restrict y,
               const real* restrict z, const real* restrict qgrid, real3 reca,
               real3 recb, real3 recc, grad_prec* restrict gx,
               grad_prec* restrict gy, grad_prec* restrict gz)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;


   real thetai1[4 * 5];
   real thetai2[4 * 5];
   real thetai3[4 * 5];
   __shared__ real sharedarray[5 * 5 * PME_BLOCKDIM];
   real* restrict array = &sharedarray[5 * 5 * threadIdx.x];


   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < n;
        ii += blockDim.x * gridDim.x) {
      real icsix = csix[ii];


      // self energy
      if CONSTEXPR (do_e) {
         int offset = ii & (bufsize - 1);
         real fs = aewald * aewald;
         fs *= fs * fs;
         fs /= 12;
         real e = fs * icsix * icsix;
         atomic_add(e, edsp, offset);
         if CONSTEXPR (do_a) {
            atomic_add(1, ndisp, offset);
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
   size_t bufsize = buffer_size();
   const auto& st = *dpme_unit;
   real aewald = st.aewald;
   int nfft1 = st.nfft1;
   int nfft2 = st.nfft2;
   int nfft3 = st.nfft3;


   assert(st.bsorder == 4);
   auto ker = edisp_cu3<Ver, 4>;
   launch_k2s(nonblk, PME_BLOCKDIM, n, ker, //
              bufsize, ndisp, edsp, csix, aewald, n, nfft1, nfft2, nfft3, x, y,
              z, st.qgrid, recipa, recipb, recipc, dedspx, dedspy, dedspz);
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
