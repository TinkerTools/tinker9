#include "add.h"
#include "epolar_trq.h"
#include "glob/spatial.h"
#include "ff/hippo/cflux.h"
#include "ff/hippo/epolar_chgpen.h"
#include "ff/image.h"
#include "launch.h"
#include "md.h"
#include "ff/pmestuf.h"
#include "seq/bsplgen.h"
#include "seq/pair_polar_chgpen.h"
#include "seq/triangle.h"
#include "ff/switch.h"
#include "tool/gpucard.h"

namespace tinker {
// ck.py Version 2.0.3
template <class Ver, class ETYP, bool CFLX>
__global__
void epolar_chgpen_cu1(int n, TINKER_IMAGE_PARAMS, count_buffer restrict np,
   energy_buffer restrict ep, virial_buffer restrict vp, grad_prec* restrict gx,
   grad_prec* restrict gy, grad_prec* restrict gz, real off, const unsigned* restrict dwinfo,
   int nexclude, const int (*restrict exclude)[2], const real (*restrict exclude_scale)[3],
   const real* restrict x, const real* restrict y, const real* restrict z,
   const Spatial::SortedAtom* restrict sorted, int nakpl, const int* restrict iakpl, int niak,
   const int* restrict iak, const int* restrict lst, real (*restrict ufld)[3],
   real (*restrict dufld)[6], const real (*restrict uind)[3], real* restrict pot,
   const real (*restrict rpole)[10], real* restrict pcore, real* restrict pval,
   const real* restrict palpha, real aewald, real f)
{
   constexpr bool do_a = Ver::a;
   constexpr bool do_e = Ver::e;
   constexpr bool do_v = Ver::v;
   constexpr bool do_g = Ver::g;
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);

   int nptl;
   if CONSTEXPR (do_a) {
      nptl = 0;
   }
   using ebuf_prec = energy_buffer_traits::type;
   ebuf_prec eptl;
   if CONSTEXPR (do_e) {
      eptl = 0;
   }
   using vbuf_prec = virial_buffer_traits::type;
   vbuf_prec vptlxx, vptlyx, vptlzx, vptlyy, vptlzy, vptlzz;
   if CONSTEXPR (do_v) {
      vptlxx = 0;
      vptlyx = 0;
      vptlzx = 0;
      vptlyy = 0;
      vptlzy = 0;
      vptlzz = 0;
   }
   __shared__ real xi[BLOCK_DIM];
   __shared__ real yi[BLOCK_DIM];
   __shared__ real zi[BLOCK_DIM];
   real xk;
   real yk;
   real zk;
   __shared__ real gxi[BLOCK_DIM];
   __shared__ real gyi[BLOCK_DIM];
   __shared__ real gzi[BLOCK_DIM];
   __shared__ real txi[BLOCK_DIM];
   __shared__ real tyi[BLOCK_DIM];
   __shared__ real tzi[BLOCK_DIM];
   __shared__ real dui0[BLOCK_DIM];
   __shared__ real dui1[BLOCK_DIM];
   __shared__ real dui2[BLOCK_DIM];
   __shared__ real dui3[BLOCK_DIM];
   __shared__ real dui4[BLOCK_DIM];
   __shared__ real dui5[BLOCK_DIM];
   __shared__ real poti[BLOCK_DIM];
   real gxk;
   real gyk;
   real gzk;
   real txk;
   real tyk;
   real tzk;
   real duk0;
   real duk1;
   real duk2;
   real duk3;
   real duk4;
   real duk5;
   real potk;
   __shared__ real ci[BLOCK_DIM];
   __shared__ real dix[BLOCK_DIM];
   __shared__ real diy[BLOCK_DIM];
   __shared__ real diz[BLOCK_DIM];
   __shared__ real qixx[BLOCK_DIM];
   __shared__ real qixy[BLOCK_DIM];
   __shared__ real qixz[BLOCK_DIM];
   __shared__ real qiyy[BLOCK_DIM];
   __shared__ real qiyz[BLOCK_DIM];
   __shared__ real qizz[BLOCK_DIM];
   __shared__ real uix[BLOCK_DIM];
   __shared__ real uiy[BLOCK_DIM];
   __shared__ real uiz[BLOCK_DIM];
   __shared__ real corei[BLOCK_DIM];
   __shared__ real alphai[BLOCK_DIM];
   __shared__ real vali[BLOCK_DIM];
   real ck;
   real dkx;
   real dky;
   real dkz;
   real qkxx;
   real qkxy;
   real qkxz;
   real qkyy;
   real qkyz;
   real qkzz;
   real ukx;
   real uky;
   real ukz;
   real corek;
   real alphak;
   real valk;

   //* /
   for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
      const int klane = threadIdx.x;
      if CONSTEXPR (do_g) {
         gxi[threadIdx.x] = 0;
         gyi[threadIdx.x] = 0;
         gzi[threadIdx.x] = 0;
         txi[threadIdx.x] = 0;
         tyi[threadIdx.x] = 0;
         tzi[threadIdx.x] = 0;
         dui0[threadIdx.x] = 0;
         dui1[threadIdx.x] = 0;
         dui2[threadIdx.x] = 0;
         dui3[threadIdx.x] = 0;
         dui4[threadIdx.x] = 0;
         dui5[threadIdx.x] = 0;
         if CONSTEXPR (CFLX)
            poti[threadIdx.x] = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
         txk = 0;
         tyk = 0;
         tzk = 0;
         duk0 = 0;
         duk1 = 0;
         duk2 = 0;
         duk3 = 0;
         duk4 = 0;
         duk5 = 0;
         if CONSTEXPR (CFLX)
            potk = 0;
      }

      int i = exclude[ii][0];
      int k = exclude[ii][1];
      real scaleb = exclude_scale[ii][1];
      real scalec = exclude_scale[ii][2];

      xi[klane] = x[i];
      yi[klane] = y[i];
      zi[klane] = z[i];
      xk = x[k];
      yk = y[k];
      zk = z[k];
      ci[klane] = rpole[i][mpl_pme_0];
      dix[klane] = rpole[i][mpl_pme_x];
      diy[klane] = rpole[i][mpl_pme_y];
      diz[klane] = rpole[i][mpl_pme_z];
      qixx[klane] = rpole[i][mpl_pme_xx];
      qixy[klane] = rpole[i][mpl_pme_xy];
      qixz[klane] = rpole[i][mpl_pme_xz];
      qiyy[klane] = rpole[i][mpl_pme_yy];
      qiyz[klane] = rpole[i][mpl_pme_yz];
      qizz[klane] = rpole[i][mpl_pme_zz];
      uix[klane] = uind[i][0];
      uiy[klane] = uind[i][1];
      uiz[klane] = uind[i][2];
      corei[klane] = pcore[i];
      alphai[klane] = palpha[i];
      vali[klane] = pval[i];
      ck = rpole[k][mpl_pme_0];
      dkx = rpole[k][mpl_pme_x];
      dky = rpole[k][mpl_pme_y];
      dkz = rpole[k][mpl_pme_z];
      qkxx = rpole[k][mpl_pme_xx];
      qkxy = rpole[k][mpl_pme_xy];
      qkxz = rpole[k][mpl_pme_xz];
      qkyy = rpole[k][mpl_pme_yy];
      qkyz = rpole[k][mpl_pme_yz];
      qkzz = rpole[k][mpl_pme_zz];
      ukx = uind[k][0];
      uky = uind[k][1];
      ukz = uind[k][2];
      corek = pcore[k];
      alphak = palpha[k];
      valk = pval[k];

      constexpr bool incl = true;
      real xr = xk - xi[klane];
      real yr = yk - yi[klane];
      real zr = zk - zi[klane];

      real e;
      real pota, potb;
      PairPolarGrad pgrad;
      zero(pgrad);

      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         pair_polar_chgpen<do_e, do_g, ETYP, CFLX>( //
            r2, xr, yr, zr, scaleb, scalec,         //
            ci[klane], dix[klane], diy[klane], diz[klane], corei[klane], vali[klane], alphai[klane],
            qixx[klane], qixy[klane], qixz[klane], qiyy[klane], qiyz[klane], qizz[klane],
            uix[klane], uiy[klane],
            uiz[klane], //
            ck, dkx, dky, dkz, corek, valk, alphak, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukx, uky,
            ukz, f, aewald, //
            e, pota, potb, pgrad);

         if CONSTEXPR (do_a)
            if (e != 0 and scaleb != 0)
               nptl += 1;
         if CONSTEXPR (do_e)
            eptl += floatTo<ebuf_prec>(e);
         if CONSTEXPR (do_g) {
            gxi[klane] += pgrad.frcx;
            gyi[klane] += pgrad.frcy;
            gzi[klane] += pgrad.frcz;
            gxk -= pgrad.frcx;
            gyk -= pgrad.frcy;
            gzk -= pgrad.frcz;

            txi[klane] += pgrad.ufldi[0];
            tyi[klane] += pgrad.ufldi[1];
            tzi[klane] += pgrad.ufldi[2];
            txk += pgrad.ufldk[0];
            tyk += pgrad.ufldk[1];
            tzk += pgrad.ufldk[2];

            dui0[klane] += pgrad.dufldi[0];
            dui1[klane] += pgrad.dufldi[1];
            dui2[klane] += pgrad.dufldi[2];
            dui3[klane] += pgrad.dufldi[3];
            dui4[klane] += pgrad.dufldi[4];
            dui5[klane] += pgrad.dufldi[5];
            duk0 += pgrad.dufldk[0];
            duk1 += pgrad.dufldk[1];
            duk2 += pgrad.dufldk[2];
            duk3 += pgrad.dufldk[3];
            duk4 += pgrad.dufldk[4];
            duk5 += pgrad.dufldk[5];
         }
         if CONSTEXPR (do_v) {
            vptlxx += floatTo<vbuf_prec>(-xr * pgrad.frcx);
            vptlyx += floatTo<vbuf_prec>(-0.5f * (yr * pgrad.frcx + xr * pgrad.frcy));
            vptlzx += floatTo<vbuf_prec>(-0.5f * (zr * pgrad.frcx + xr * pgrad.frcz));
            vptlyy += floatTo<vbuf_prec>(-yr * pgrad.frcy);
            vptlzy += floatTo<vbuf_prec>(-0.5f * (zr * pgrad.frcy + yr * pgrad.frcz));
            vptlzz += floatTo<vbuf_prec>(-zr * pgrad.frcz);
         }
         if CONSTEXPR (CFLX) {
            poti[klane] += pota;
            potk += potb;
         }
      } // end if (include)

      if CONSTEXPR (do_g) {
         atomic_add(gxi[threadIdx.x], gx, i);
         atomic_add(gyi[threadIdx.x], gy, i);
         atomic_add(gzi[threadIdx.x], gz, i);
         atomic_add(txi[threadIdx.x], &ufld[i][0]);
         atomic_add(tyi[threadIdx.x], &ufld[i][1]);
         atomic_add(tzi[threadIdx.x], &ufld[i][2]);
         atomic_add(dui0[threadIdx.x], &dufld[i][0]);
         atomic_add(dui1[threadIdx.x], &dufld[i][1]);
         atomic_add(dui2[threadIdx.x], &dufld[i][2]);
         atomic_add(dui3[threadIdx.x], &dufld[i][3]);
         atomic_add(dui4[threadIdx.x], &dufld[i][4]);
         atomic_add(dui5[threadIdx.x], &dufld[i][5]);
         if CONSTEXPR (CFLX)
            atomic_add(poti[threadIdx.x], pot, i);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
         atomic_add(txk, &ufld[k][0]);
         atomic_add(tyk, &ufld[k][1]);
         atomic_add(tzk, &ufld[k][2]);
         atomic_add(duk0, &dufld[k][0]);
         atomic_add(duk1, &dufld[k][1]);
         atomic_add(duk2, &dufld[k][2]);
         atomic_add(duk3, &dufld[k][3]);
         atomic_add(duk4, &dufld[k][4]);
         atomic_add(duk5, &dufld[k][5]);
         if CONSTEXPR (CFLX)
            atomic_add(potk, pot, k);
      }
   }
   // */

   for (int iw = iwarp; iw < nakpl; iw += nwarp) {
      if CONSTEXPR (do_g) {
         gxi[threadIdx.x] = 0;
         gyi[threadIdx.x] = 0;
         gzi[threadIdx.x] = 0;
         txi[threadIdx.x] = 0;
         tyi[threadIdx.x] = 0;
         tzi[threadIdx.x] = 0;
         dui0[threadIdx.x] = 0;
         dui1[threadIdx.x] = 0;
         dui2[threadIdx.x] = 0;
         dui3[threadIdx.x] = 0;
         dui4[threadIdx.x] = 0;
         dui5[threadIdx.x] = 0;
         if CONSTEXPR (CFLX)
            poti[threadIdx.x] = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
         txk = 0;
         tyk = 0;
         tzk = 0;
         duk0 = 0;
         duk1 = 0;
         duk2 = 0;
         duk3 = 0;
         duk4 = 0;
         duk5 = 0;
         if CONSTEXPR (CFLX)
            potk = 0;
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
      xi[threadIdx.x] = sorted[atomi].x;
      yi[threadIdx.x] = sorted[atomi].y;
      zi[threadIdx.x] = sorted[atomi].z;
      xk = sorted[atomk].x;
      yk = sorted[atomk].y;
      zk = sorted[atomk].z;

      ci[threadIdx.x] = rpole[i][mpl_pme_0];
      dix[threadIdx.x] = rpole[i][mpl_pme_x];
      diy[threadIdx.x] = rpole[i][mpl_pme_y];
      diz[threadIdx.x] = rpole[i][mpl_pme_z];
      qixx[threadIdx.x] = rpole[i][mpl_pme_xx];
      qixy[threadIdx.x] = rpole[i][mpl_pme_xy];
      qixz[threadIdx.x] = rpole[i][mpl_pme_xz];
      qiyy[threadIdx.x] = rpole[i][mpl_pme_yy];
      qiyz[threadIdx.x] = rpole[i][mpl_pme_yz];
      qizz[threadIdx.x] = rpole[i][mpl_pme_zz];
      uix[threadIdx.x] = uind[i][0];
      uiy[threadIdx.x] = uind[i][1];
      uiz[threadIdx.x] = uind[i][2];
      corei[threadIdx.x] = pcore[i];
      alphai[threadIdx.x] = palpha[i];
      vali[threadIdx.x] = pval[i];
      ck = rpole[k][mpl_pme_0];
      dkx = rpole[k][mpl_pme_x];
      dky = rpole[k][mpl_pme_y];
      dkz = rpole[k][mpl_pme_z];
      qkxx = rpole[k][mpl_pme_xx];
      qkxy = rpole[k][mpl_pme_xy];
      qkxz = rpole[k][mpl_pme_xz];
      qkyy = rpole[k][mpl_pme_yy];
      qkyz = rpole[k][mpl_pme_yz];
      qkzz = rpole[k][mpl_pme_zz];
      ukx = uind[k][0];
      uky = uind[k][1];
      ukz = uind[k][2];
      corek = pcore[k];
      alphak = palpha[k];
      valk = pval[k];

      unsigned int dwinfo0 = dwinfo[iw * WARP_SIZE + ilane];
      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         bool incl = iid < kid and kid < n;
         int srcmask = 1 << srclane;
         incl = incl and (dwinfo0 & srcmask) == 0;
         real xr = xk - xi[klane];
         real yr = yk - yi[klane];
         real zr = zk - zi[klane];

         real e;
         real pota, potb;
         PairPolarGrad pgrad;
         zero(pgrad);

         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_polar_chgpen<do_e, do_g, ETYP, CFLX>( //
               r2, xr, yr, zr, 1, 1,                   //
               ci[klane], dix[klane], diy[klane], diz[klane], corei[klane], vali[klane],
               alphai[klane], qixx[klane], qixy[klane], qixz[klane], qiyy[klane], qiyz[klane],
               qizz[klane], uix[klane], uiy[klane], uiz[klane], //
               ck, dkx, dky, dkz, corek, valk, alphak, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukx, uky,
               ukz, f, aewald, //
               e, pota, potb, pgrad);

            if CONSTEXPR (do_a)
               if (e != 0)
                  nptl += 1;
            if CONSTEXPR (do_e)
               eptl += floatTo<ebuf_prec>(e);
            if CONSTEXPR (do_g) {
               gxi[klane] += pgrad.frcx;
               gyi[klane] += pgrad.frcy;
               gzi[klane] += pgrad.frcz;
               gxk -= pgrad.frcx;
               gyk -= pgrad.frcy;
               gzk -= pgrad.frcz;

               txi[klane] += pgrad.ufldi[0];
               tyi[klane] += pgrad.ufldi[1];
               tzi[klane] += pgrad.ufldi[2];
               txk += pgrad.ufldk[0];
               tyk += pgrad.ufldk[1];
               tzk += pgrad.ufldk[2];

               dui0[klane] += pgrad.dufldi[0];
               dui1[klane] += pgrad.dufldi[1];
               dui2[klane] += pgrad.dufldi[2];
               dui3[klane] += pgrad.dufldi[3];
               dui4[klane] += pgrad.dufldi[4];
               dui5[klane] += pgrad.dufldi[5];
               duk0 += pgrad.dufldk[0];
               duk1 += pgrad.dufldk[1];
               duk2 += pgrad.dufldk[2];
               duk3 += pgrad.dufldk[3];
               duk4 += pgrad.dufldk[4];
               duk5 += pgrad.dufldk[5];
            }
            if CONSTEXPR (do_v) {
               vptlxx += floatTo<vbuf_prec>(-xr * pgrad.frcx);
               vptlyx += floatTo<vbuf_prec>(-0.5f * (yr * pgrad.frcx + xr * pgrad.frcy));
               vptlzx += floatTo<vbuf_prec>(-0.5f * (zr * pgrad.frcx + xr * pgrad.frcz));
               vptlyy += floatTo<vbuf_prec>(-yr * pgrad.frcy);
               vptlzy += floatTo<vbuf_prec>(-0.5f * (zr * pgrad.frcy + yr * pgrad.frcz));
               vptlzz += floatTo<vbuf_prec>(-zr * pgrad.frcz);
            }
            if CONSTEXPR (CFLX) {
               poti[klane] += pota;
               potk += potb;
            }
         } // end if (include)

         iid = __shfl_sync(ALL_LANES, iid, ilane + 1);
      }

      if CONSTEXPR (do_g) {
         atomic_add(gxi[threadIdx.x], gx, i);
         atomic_add(gyi[threadIdx.x], gy, i);
         atomic_add(gzi[threadIdx.x], gz, i);
         atomic_add(txi[threadIdx.x], &ufld[i][0]);
         atomic_add(tyi[threadIdx.x], &ufld[i][1]);
         atomic_add(tzi[threadIdx.x], &ufld[i][2]);
         atomic_add(dui0[threadIdx.x], &dufld[i][0]);
         atomic_add(dui1[threadIdx.x], &dufld[i][1]);
         atomic_add(dui2[threadIdx.x], &dufld[i][2]);
         atomic_add(dui3[threadIdx.x], &dufld[i][3]);
         atomic_add(dui4[threadIdx.x], &dufld[i][4]);
         atomic_add(dui5[threadIdx.x], &dufld[i][5]);
         if CONSTEXPR (CFLX)
            atomic_add(poti[threadIdx.x], pot, i);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
         atomic_add(txk, &ufld[k][0]);
         atomic_add(tyk, &ufld[k][1]);
         atomic_add(tzk, &ufld[k][2]);
         atomic_add(duk0, &dufld[k][0]);
         atomic_add(duk1, &dufld[k][1]);
         atomic_add(duk2, &dufld[k][2]);
         atomic_add(duk3, &dufld[k][3]);
         atomic_add(duk4, &dufld[k][4]);
         atomic_add(duk5, &dufld[k][5]);
         if CONSTEXPR (CFLX)
            atomic_add(potk, pot, k);
      }
   }

   for (int iw = iwarp; iw < niak; iw += nwarp) {
      if CONSTEXPR (do_g) {
         gxi[threadIdx.x] = 0;
         gyi[threadIdx.x] = 0;
         gzi[threadIdx.x] = 0;
         txi[threadIdx.x] = 0;
         tyi[threadIdx.x] = 0;
         tzi[threadIdx.x] = 0;
         dui0[threadIdx.x] = 0;
         dui1[threadIdx.x] = 0;
         dui2[threadIdx.x] = 0;
         dui3[threadIdx.x] = 0;
         dui4[threadIdx.x] = 0;
         dui5[threadIdx.x] = 0;
         if CONSTEXPR (CFLX)
            poti[threadIdx.x] = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
         txk = 0;
         tyk = 0;
         tzk = 0;
         duk0 = 0;
         duk1 = 0;
         duk2 = 0;
         duk3 = 0;
         duk4 = 0;
         duk5 = 0;
         if CONSTEXPR (CFLX)
            potk = 0;
      }

      int ty = iak[iw];
      int atomi = ty * WARP_SIZE + ilane;
      int i = sorted[atomi].unsorted;
      int atomk = lst[iw * WARP_SIZE + ilane];
      int k = sorted[atomk].unsorted;
      xi[threadIdx.x] = sorted[atomi].x;
      yi[threadIdx.x] = sorted[atomi].y;
      zi[threadIdx.x] = sorted[atomi].z;
      xk = sorted[atomk].x;
      yk = sorted[atomk].y;
      zk = sorted[atomk].z;

      ci[threadIdx.x] = rpole[i][mpl_pme_0];
      dix[threadIdx.x] = rpole[i][mpl_pme_x];
      diy[threadIdx.x] = rpole[i][mpl_pme_y];
      diz[threadIdx.x] = rpole[i][mpl_pme_z];
      qixx[threadIdx.x] = rpole[i][mpl_pme_xx];
      qixy[threadIdx.x] = rpole[i][mpl_pme_xy];
      qixz[threadIdx.x] = rpole[i][mpl_pme_xz];
      qiyy[threadIdx.x] = rpole[i][mpl_pme_yy];
      qiyz[threadIdx.x] = rpole[i][mpl_pme_yz];
      qizz[threadIdx.x] = rpole[i][mpl_pme_zz];
      uix[threadIdx.x] = uind[i][0];
      uiy[threadIdx.x] = uind[i][1];
      uiz[threadIdx.x] = uind[i][2];
      corei[threadIdx.x] = pcore[i];
      alphai[threadIdx.x] = palpha[i];
      vali[threadIdx.x] = pval[i];
      ck = rpole[k][mpl_pme_0];
      dkx = rpole[k][mpl_pme_x];
      dky = rpole[k][mpl_pme_y];
      dkz = rpole[k][mpl_pme_z];
      qkxx = rpole[k][mpl_pme_xx];
      qkxy = rpole[k][mpl_pme_xy];
      qkxz = rpole[k][mpl_pme_xz];
      qkyy = rpole[k][mpl_pme_yy];
      qkyz = rpole[k][mpl_pme_yz];
      qkzz = rpole[k][mpl_pme_zz];
      ukx = uind[k][0];
      uky = uind[k][1];
      ukz = uind[k][2];
      corek = pcore[k];
      alphak = palpha[k];
      valk = pval[k];

      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         bool incl = atomk > 0;
         real xr = xk - xi[klane];
         real yr = yk - yi[klane];
         real zr = zk - zi[klane];

         real e;
         real pota, potb;
         PairPolarGrad pgrad;
         zero(pgrad);

         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_polar_chgpen<do_e, do_g, ETYP, CFLX>( //
               r2, xr, yr, zr, 1, 1,                   //
               ci[klane], dix[klane], diy[klane], diz[klane], corei[klane], vali[klane],
               alphai[klane], qixx[klane], qixy[klane], qixz[klane], qiyy[klane], qiyz[klane],
               qizz[klane], uix[klane], uiy[klane], uiz[klane], //
               ck, dkx, dky, dkz, corek, valk, alphak, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukx, uky,
               ukz, f, aewald, //
               e, pota, potb, pgrad);

            if CONSTEXPR (do_a)
               if (e != 0)
                  nptl += 1;
            if CONSTEXPR (do_e)
               eptl += floatTo<ebuf_prec>(e);
            if CONSTEXPR (do_g) {
               gxi[klane] += pgrad.frcx;
               gyi[klane] += pgrad.frcy;
               gzi[klane] += pgrad.frcz;
               gxk -= pgrad.frcx;
               gyk -= pgrad.frcy;
               gzk -= pgrad.frcz;

               txi[klane] += pgrad.ufldi[0];
               tyi[klane] += pgrad.ufldi[1];
               tzi[klane] += pgrad.ufldi[2];
               txk += pgrad.ufldk[0];
               tyk += pgrad.ufldk[1];
               tzk += pgrad.ufldk[2];

               dui0[klane] += pgrad.dufldi[0];
               dui1[klane] += pgrad.dufldi[1];
               dui2[klane] += pgrad.dufldi[2];
               dui3[klane] += pgrad.dufldi[3];
               dui4[klane] += pgrad.dufldi[4];
               dui5[klane] += pgrad.dufldi[5];
               duk0 += pgrad.dufldk[0];
               duk1 += pgrad.dufldk[1];
               duk2 += pgrad.dufldk[2];
               duk3 += pgrad.dufldk[3];
               duk4 += pgrad.dufldk[4];
               duk5 += pgrad.dufldk[5];
            }
            if CONSTEXPR (do_v) {
               vptlxx += floatTo<vbuf_prec>(-xr * pgrad.frcx);
               vptlyx += floatTo<vbuf_prec>(-0.5f * (yr * pgrad.frcx + xr * pgrad.frcy));
               vptlzx += floatTo<vbuf_prec>(-0.5f * (zr * pgrad.frcx + xr * pgrad.frcz));
               vptlyy += floatTo<vbuf_prec>(-yr * pgrad.frcy);
               vptlzy += floatTo<vbuf_prec>(-0.5f * (zr * pgrad.frcy + yr * pgrad.frcz));
               vptlzz += floatTo<vbuf_prec>(-zr * pgrad.frcz);
            }
            if CONSTEXPR (CFLX) {
               poti[klane] += pota;
               potk += potb;
            }
         } // end if (include)
      }

      if CONSTEXPR (do_g) {
         atomic_add(gxi[threadIdx.x], gx, i);
         atomic_add(gyi[threadIdx.x], gy, i);
         atomic_add(gzi[threadIdx.x], gz, i);
         atomic_add(txi[threadIdx.x], &ufld[i][0]);
         atomic_add(tyi[threadIdx.x], &ufld[i][1]);
         atomic_add(tzi[threadIdx.x], &ufld[i][2]);
         atomic_add(dui0[threadIdx.x], &dufld[i][0]);
         atomic_add(dui1[threadIdx.x], &dufld[i][1]);
         atomic_add(dui2[threadIdx.x], &dufld[i][2]);
         atomic_add(dui3[threadIdx.x], &dufld[i][3]);
         atomic_add(dui4[threadIdx.x], &dufld[i][4]);
         atomic_add(dui5[threadIdx.x], &dufld[i][5]);
         if CONSTEXPR (CFLX)
            atomic_add(poti[threadIdx.x], pot, i);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
         atomic_add(txk, &ufld[k][0]);
         atomic_add(tyk, &ufld[k][1]);
         atomic_add(tzk, &ufld[k][2]);
         atomic_add(duk0, &dufld[k][0]);
         atomic_add(duk1, &dufld[k][1]);
         atomic_add(duk2, &dufld[k][2]);
         atomic_add(duk3, &dufld[k][3]);
         atomic_add(duk4, &dufld[k][4]);
         atomic_add(duk5, &dufld[k][5]);
         if CONSTEXPR (CFLX)
            atomic_add(potk, pot, k);
      }
   }

   if CONSTEXPR (do_a) {
      atomic_add(nptl, np, ithread);
   }
   if CONSTEXPR (do_e) {
      atomic_add(eptl, ep, ithread);
   }
   if CONSTEXPR (do_v) {
      atomic_add(vptlxx, vptlyx, vptlzx, vptlyy, vptlzy, vptlzz, vp, ithread);
   }
}

template <class Ver, class ETYP, int CFLX>
void epolar_chgpen_cu(const real (*uind)[3])
{
   constexpr bool do_g = Ver::g;
   const auto& st = *mspatial_v2_unit;
   real off;

   if CONSTEXPR (eq<ETYP, EWALD>())
      off = switch_off(switch_ewald);
   else
      off = switch_off(switch_mpole);

   const real f = 0.5f * electric / dielec;
   real aewald = 0;
   if CONSTEXPR (eq<ETYP, EWALD>()) {
      PMEUnit pu = ppme_unit;
      aewald = pu->aewald;
   }

   if CONSTEXPR (do_g)
      darray::zero(g::q0, n, ufld, dufld);

   int ngrid = gpuGridSize(BLOCK_DIM);
   epolar_chgpen_cu1<Ver, ETYP, CFLX><<<ngrid, BLOCK_DIM, 0, g::s0>>>(st.n, TINKER_IMAGE_ARGS, nep,
      ep, vir_ep, depx, depy, depz, off, st.si1.bit0, nmdwexclude, mdwexclude, mdwexclude_scale,
      st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst, ufld, dufld, uind,
      pot, rpole, pcore, pval, palpha, aewald, f);

   // torque
   if CONSTEXPR (do_g) {
      launch_k1s(g::s0, n, epolar_trq_cu, //
         trqx, trqy, trqz, n, rpole, ufld, dufld);
   }
}

void epolar_chgpen_nonewald_cu(int vers, int use_cf, const real (*uind)[3])
{
   if (use_cf) {
      if (vers == calc::v0) {
         // epolar_chgpen_cu<calc::V0, NON_EWALD, 1>(uind);
         assert(false && "CFLX must compute gradient.");
      } else if (vers == calc::v1) {
         epolar_chgpen_cu<calc::V1, NON_EWALD, 1>(uind);
      } else if (vers == calc::v3) {
         // epolar_chgpen_cu<calc::V3, NON_EWALD, 1>(uind);
         assert(false && "CFLX must compute gradient.");
      } else if (vers == calc::v4) {
         epolar_chgpen_cu<calc::V4, NON_EWALD, 1>(uind);
      } else if (vers == calc::v5) {
         epolar_chgpen_cu<calc::V5, NON_EWALD, 1>(uind);
      } else if (vers == calc::v6) {
         epolar_chgpen_cu<calc::V6, NON_EWALD, 1>(uind);
      }
   } else {
      if (vers == calc::v0) {
         epolar_chgpen_cu<calc::V0, NON_EWALD, 0>(uind);
      } else if (vers == calc::v1) {
         epolar_chgpen_cu<calc::V1, NON_EWALD, 0>(uind);
      } else if (vers == calc::v3) {
         epolar_chgpen_cu<calc::V3, NON_EWALD, 0>(uind);
      } else if (vers == calc::v4) {
         epolar_chgpen_cu<calc::V4, NON_EWALD, 0>(uind);
      } else if (vers == calc::v5) {
         epolar_chgpen_cu<calc::V5, NON_EWALD, 0>(uind);
      } else if (vers == calc::v6) {
         epolar_chgpen_cu<calc::V6, NON_EWALD, 0>(uind);
      }
   }
}

void epolar_chgpen_ewald_real_cu(int vers, int use_cf, const real (*uind)[3])
{
   if (use_cf) {
      if (vers == calc::v0) {
         // epolar_chgpen_cu<calc::V0, EWALD, 1>(uind);
         assert(false && "CFLX must compute gradient.");
      } else if (vers == calc::v1) {
         epolar_chgpen_cu<calc::V1, EWALD, 1>(uind);
      } else if (vers == calc::v3) {
         // epolar_chgpen_cu<calc::V3, EWALD, 1>(uind);
         assert(false && "CFLX must compute gradient.");
      } else if (vers == calc::v4) {
         epolar_chgpen_cu<calc::V4, EWALD, 1>(uind);
      } else if (vers == calc::v5) {
         epolar_chgpen_cu<calc::V5, EWALD, 1>(uind);
      } else if (vers == calc::v6) {
         epolar_chgpen_cu<calc::V6, EWALD, 1>(uind);
      }
   } else {
      if (vers == calc::v0) {
         epolar_chgpen_cu<calc::V0, EWALD, 0>(uind);
      } else if (vers == calc::v1) {
         epolar_chgpen_cu<calc::V1, EWALD, 0>(uind);
      } else if (vers == calc::v3) {
         epolar_chgpen_cu<calc::V3, EWALD, 0>(uind);
      } else if (vers == calc::v4) {
         epolar_chgpen_cu<calc::V4, EWALD, 0>(uind);
      } else if (vers == calc::v5) {
         epolar_chgpen_cu<calc::V5, EWALD, 0>(uind);
      } else if (vers == calc::v6) {
         epolar_chgpen_cu<calc::V6, EWALD, 0>(uind);
      }
   }
}
}
