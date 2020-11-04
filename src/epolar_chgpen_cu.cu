#include "add.h"
#include "cflux.h"
#include "epolar_chgpen.h"
#include "epolar_trq.h"
#include "glob.spatial.h"
#include "image.h"
#include "launch.h"
#include "md.h"
#include "pmestuf.h"
#include "seq_bsplgen.h"
#include "seq_pair_polar_chgpen.h"
#include "seq_triangle.h"
#include "switch.h"
#include "tool/gpu_card.h"


namespace tinker {
template <class Ver, class ETYP, int CFLX>
__global__
void epolar_chgpen_cu1(
   int n, TINKER_IMAGE_PARAMS, count_buffer restrict np,
   energy_buffer restrict ep, virial_buffer restrict vp, grad_prec* restrict gx,
   grad_prec* restrict gy, grad_prec* restrict gz, real off,
   const unsigned* restrict dwinfo, int nexclude,
   const int (*restrict exclude)[2], const real (*restrict exclude_scale)[3],
   const real* restrict x, const real* restrict y, const real* restrict z,
   const Spatial::SortedAtom* restrict sorted, int nakpl,
   const int* restrict iakpl, int niak, const int* restrict iak,
   const int* restrict lst, real (*restrict ufld)[3], real (*restrict dufld)[6],
   const real (*restrict uind)[3], real* restrict pot,
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


   __shared__ real shxi[BLOCK_DIM];
   __shared__ real shyi[BLOCK_DIM];
   __shared__ real shzi[BLOCK_DIM];
   real xk;
   real yk;
   real zk;
   __shared__ real shgxi[BLOCK_DIM];
   __shared__ real shgyi[BLOCK_DIM];
   __shared__ real shgzi[BLOCK_DIM];
   __shared__ real shtxi[BLOCK_DIM];
   __shared__ real shtyi[BLOCK_DIM];
   __shared__ real shtzi[BLOCK_DIM];
   __shared__ real shdui0[BLOCK_DIM];
   __shared__ real shdui1[BLOCK_DIM];
   __shared__ real shdui2[BLOCK_DIM];
   __shared__ real shdui3[BLOCK_DIM];
   __shared__ real shdui4[BLOCK_DIM];
   __shared__ real shdui5[BLOCK_DIM];
   __shared__ real shpoti[BLOCK_DIM];
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
   __shared__ real shci[BLOCK_DIM];
   __shared__ real shdix[BLOCK_DIM];
   __shared__ real shdiy[BLOCK_DIM];
   __shared__ real shdiz[BLOCK_DIM];
   __shared__ real shqixx[BLOCK_DIM];
   __shared__ real shqixy[BLOCK_DIM];
   __shared__ real shqixz[BLOCK_DIM];
   __shared__ real shqiyy[BLOCK_DIM];
   __shared__ real shqiyz[BLOCK_DIM];
   __shared__ real shqizz[BLOCK_DIM];
   __shared__ real shuix[BLOCK_DIM];
   __shared__ real shuiy[BLOCK_DIM];
   __shared__ real shuiz[BLOCK_DIM];
   __shared__ real shcorei[BLOCK_DIM];
   __shared__ real shalphai[BLOCK_DIM];
   __shared__ real shvali[BLOCK_DIM];
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
   // exclude
   for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
      if CONSTEXPR (do_g) {
         shgxi[threadIdx.x] = 0;
         shgyi[threadIdx.x] = 0;
         shgzi[threadIdx.x] = 0;
         shtxi[threadIdx.x] = 0;
         shtyi[threadIdx.x] = 0;
         shtzi[threadIdx.x] = 0;
         shdui0[threadIdx.x] = 0;
         shdui1[threadIdx.x] = 0;
         shdui2[threadIdx.x] = 0;
         shdui3[threadIdx.x] = 0;
         shdui4[threadIdx.x] = 0;
         shdui5[threadIdx.x] = 0;
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
      }
      if CONSTEXPR (CFLX) {
         shpoti[threadIdx.x] = 0;
         potk = 0;
      }


      int shi = exclude[ii][0];
      int k = exclude[ii][1];
      real scaleb = exclude_scale[ii][1];
      real scalec = exclude_scale[ii][2];


      real xi = x[shi];
      real yi = y[shi];
      real zi = z[shi];
      xk = x[k];
      yk = y[k];
      zk = z[k];
      real ci = rpole[shi][mpl_pme_0];
      real dix = rpole[shi][mpl_pme_x];
      real diy = rpole[shi][mpl_pme_y];
      real diz = rpole[shi][mpl_pme_z];
      real qixx = rpole[shi][mpl_pme_xx];
      real qixy = rpole[shi][mpl_pme_xy];
      real qixz = rpole[shi][mpl_pme_xz];
      real qiyy = rpole[shi][mpl_pme_yy];
      real qiyz = rpole[shi][mpl_pme_yz];
      real qizz = rpole[shi][mpl_pme_zz];
      real uix = uind[shi][0];
      real uiy = uind[shi][1];
      real uiz = uind[shi][2];
      real corei = pcore[shi];
      real alphai = palpha[shi];
      real vali = pval[shi];
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


      int klane = threadIdx.x;
      constexpr bool incl = true;
      real xr = xk - xi;
      real yr = yk - yi;
      real zr = zk - zi;

      real e = 0;
      real pota, potb;
      PairPolarGrad pgrad;
      zero(pgrad);

      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         pair_polar_chgpen<do_e, do_g, ETYP, CFLX>( //
            r2, xr, yr, zr, scaleb, scalec,         //
            ci, dix, diy, diz, corei, vali, alphai, qixx, qixy, qixz, qiyy,
            qiyz, qizz, uix, uiy, uiz, //
            ck, dkx, dky, dkz, corek, valk, alphak, qkxx, qkxy, qkxz, qkyy,
            qkyz, qkzz, ukx, uky, ukz, f, aewald, //
            e, pota, potb, pgrad);

         if CONSTEXPR (do_a)
            if (e != 0 and scaleb != 0)
               nptl += 1;
         if CONSTEXPR (do_e)
            eptl += cvt_to<ebuf_prec>(e);
         if CONSTEXPR (do_g) {
            shgxi[klane] += pgrad.frcx;
            shgyi[klane] += pgrad.frcy;
            shgzi[klane] += pgrad.frcz;
            gxk -= pgrad.frcx;
            gyk -= pgrad.frcy;
            gzk -= pgrad.frcz;


            shtxi[klane] += pgrad.ufldi[0];
            shtyi[klane] += pgrad.ufldi[1];
            shtzi[klane] += pgrad.ufldi[2];
            txk += pgrad.ufldk[0];
            tyk += pgrad.ufldk[1];
            tzk += pgrad.ufldk[2];


            shdui0[klane] += pgrad.dufldi[0];
            shdui1[klane] += pgrad.dufldi[1];
            shdui2[klane] += pgrad.dufldi[2];
            shdui3[klane] += pgrad.dufldi[3];
            shdui4[klane] += pgrad.dufldi[4];
            shdui5[klane] += pgrad.dufldi[5];
            duk0 += pgrad.dufldk[0];
            duk1 += pgrad.dufldk[1];
            duk2 += pgrad.dufldk[2];
            duk3 += pgrad.dufldk[3];
            duk4 += pgrad.dufldk[4];
            duk5 += pgrad.dufldk[5];
         }
         if CONSTEXPR (do_v) {
            vptlxx += cvt_to<vbuf_prec>(-xr * pgrad.frcx);
            vptlyx +=
               cvt_to<vbuf_prec>(-0.5f * (yr * pgrad.frcx + xr * pgrad.frcy));
            vptlzx +=
               cvt_to<vbuf_prec>(-0.5f * (zr * pgrad.frcx + xr * pgrad.frcz));
            vptlyy += cvt_to<vbuf_prec>(-yr * pgrad.frcy);
            vptlzy +=
               cvt_to<vbuf_prec>(-0.5f * (zr * pgrad.frcy + yr * pgrad.frcz));
            vptlzz += cvt_to<vbuf_prec>(-zr * pgrad.frcz);
         }
         if CONSTEXPR (CFLX) {
            shpoti[klane] += pota;
            potk += potb;
         }
      } // end if (include)


      if CONSTEXPR (do_g) {
         atomic_add(shgxi[threadIdx.x], gx, shi);
         atomic_add(shgyi[threadIdx.x], gy, shi);
         atomic_add(shgzi[threadIdx.x], gz, shi);
         atomic_add(shtxi[threadIdx.x], &ufld[shi][0]);
         atomic_add(shtyi[threadIdx.x], &ufld[shi][1]);
         atomic_add(shtzi[threadIdx.x], &ufld[shi][2]);
         atomic_add(shdui0[threadIdx.x], &dufld[shi][0]);
         atomic_add(shdui1[threadIdx.x], &dufld[shi][1]);
         atomic_add(shdui2[threadIdx.x], &dufld[shi][2]);
         atomic_add(shdui3[threadIdx.x], &dufld[shi][3]);
         atomic_add(shdui4[threadIdx.x], &dufld[shi][4]);
         atomic_add(shdui5[threadIdx.x], &dufld[shi][5]);
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
      }
      if CONSTEXPR (CFLX) {
         atomic_add(shpoti[threadIdx.x], pot, shi);
         atomic_add(potk, pot, k);
      }
   }
   // */


   //* /
   // block pairs that have scale factors
   for (int iw = iwarp; iw < nakpl; iw += nwarp) {
      if CONSTEXPR (do_g) {
         shgxi[threadIdx.x] = 0;
         shgyi[threadIdx.x] = 0;
         shgzi[threadIdx.x] = 0;
         shtxi[threadIdx.x] = 0;
         shtyi[threadIdx.x] = 0;
         shtzi[threadIdx.x] = 0;
         shdui0[threadIdx.x] = 0;
         shdui1[threadIdx.x] = 0;
         shdui2[threadIdx.x] = 0;
         shdui3[threadIdx.x] = 0;
         shdui4[threadIdx.x] = 0;
         shdui5[threadIdx.x] = 0;
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
      }
      if CONSTEXPR (CFLX) {
         shpoti[threadIdx.x] = 0;
         potk = 0;
      }


      int tri, tx, ty;
      tri = iakpl[iw];
      tri_to_xy(tri, tx, ty);


      int shiid = ty * WARP_SIZE + ilane;
      int shatomi = min(shiid, n - 1);
      int shi = sorted[shatomi].unsorted;
      int kid = tx * WARP_SIZE + ilane;
      int atomk = min(kid, n - 1);
      int k = sorted[atomk].unsorted;
      shxi[threadIdx.x] = sorted[shatomi].x;
      shyi[threadIdx.x] = sorted[shatomi].y;
      shzi[threadIdx.x] = sorted[shatomi].z;
      xk = sorted[atomk].x;
      yk = sorted[atomk].y;
      zk = sorted[atomk].z;


      shci[threadIdx.x] = rpole[shi][mpl_pme_0];
      shdix[threadIdx.x] = rpole[shi][mpl_pme_x];
      shdiy[threadIdx.x] = rpole[shi][mpl_pme_y];
      shdiz[threadIdx.x] = rpole[shi][mpl_pme_z];
      shqixx[threadIdx.x] = rpole[shi][mpl_pme_xx];
      shqixy[threadIdx.x] = rpole[shi][mpl_pme_xy];
      shqixz[threadIdx.x] = rpole[shi][mpl_pme_xz];
      shqiyy[threadIdx.x] = rpole[shi][mpl_pme_yy];
      shqiyz[threadIdx.x] = rpole[shi][mpl_pme_yz];
      shqizz[threadIdx.x] = rpole[shi][mpl_pme_zz];
      shuix[threadIdx.x] = uind[shi][0];
      shuiy[threadIdx.x] = uind[shi][1];
      shuiz[threadIdx.x] = uind[shi][2];
      shcorei[threadIdx.x] = pcore[shi];
      shalphai[threadIdx.x] = palpha[shi];
      shvali[threadIdx.x] = pval[shi];
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
         int srcmask = 1 << srclane;
         int klane = srclane + threadIdx.x - ilane;
         int iid = shiid;
         real xi = shxi[klane];
         real yi = shyi[klane];
         real zi = shzi[klane];
         real ci = shci[klane];
         real dix = shdix[klane];
         real diy = shdiy[klane];
         real diz = shdiz[klane];
         real qixx = shqixx[klane];
         real qixy = shqixy[klane];
         real qixz = shqixz[klane];
         real qiyy = shqiyy[klane];
         real qiyz = shqiyz[klane];
         real qizz = shqizz[klane];
         real uix = shuix[klane];
         real uiy = shuiy[klane];
         real uiz = shuiz[klane];
         real corei = shcorei[klane];
         real alphai = shalphai[klane];
         real vali = shvali[klane];


         bool incl = iid < kid and kid < n;
         incl = incl and (dwinfo0 & srcmask) == 0;
         real scaleb = 1;
         real scalec = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;

         real e;
         real pota = 0, potb = 0;
         PairPolarGrad pgrad;
         zero(pgrad);

         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_polar_chgpen<do_e, do_g, ETYP, CFLX>( //
               r2, xr, yr, zr, scaleb, scalec,         //
               ci, dix, diy, diz, corei, vali, alphai, qixx, qixy, qixz, qiyy,
               qiyz, qizz, uix, uiy, uiz, //
               ck, dkx, dky, dkz, corek, valk, alphak, qkxx, qkxy, qkxz, qkyy,
               qkyz, qkzz, ukx, uky, ukz, f, aewald, //
               e, pota, potb, pgrad);

            if CONSTEXPR (do_a)
               if (e != 0)
                  nptl += 1;
            if CONSTEXPR (do_e)
               eptl += cvt_to<ebuf_prec>(e);
            if CONSTEXPR (do_g) {
               shgxi[klane] += pgrad.frcx;
               shgyi[klane] += pgrad.frcy;
               shgzi[klane] += pgrad.frcz;
               gxk -= pgrad.frcx;
               gyk -= pgrad.frcy;
               gzk -= pgrad.frcz;


               shtxi[klane] += pgrad.ufldi[0];
               shtyi[klane] += pgrad.ufldi[1];
               shtzi[klane] += pgrad.ufldi[2];
               txk += pgrad.ufldk[0];
               tyk += pgrad.ufldk[1];
               tzk += pgrad.ufldk[2];


               shdui0[klane] += pgrad.dufldi[0];
               shdui1[klane] += pgrad.dufldi[1];
               shdui2[klane] += pgrad.dufldi[2];
               shdui3[klane] += pgrad.dufldi[3];
               shdui4[klane] += pgrad.dufldi[4];
               shdui5[klane] += pgrad.dufldi[5];
               duk0 += pgrad.dufldk[0];
               duk1 += pgrad.dufldk[1];
               duk2 += pgrad.dufldk[2];
               duk3 += pgrad.dufldk[3];
               duk4 += pgrad.dufldk[4];
               duk5 += pgrad.dufldk[5];
            }
            if CONSTEXPR (do_v) {
               vptlxx += cvt_to<vbuf_prec>(-xr * pgrad.frcx);
               vptlyx += cvt_to<vbuf_prec>(-0.5f *
                                           (yr * pgrad.frcx + xr * pgrad.frcy));
               vptlzx += cvt_to<vbuf_prec>(-0.5f *
                                           (zr * pgrad.frcx + xr * pgrad.frcz));
               vptlyy += cvt_to<vbuf_prec>(-yr * pgrad.frcy);
               vptlzy += cvt_to<vbuf_prec>(-0.5f *
                                           (zr * pgrad.frcy + yr * pgrad.frcz));
               vptlzz += cvt_to<vbuf_prec>(-zr * pgrad.frcz);
            }
            if CONSTEXPR (CFLX) {
               shpoti[klane] += pota;
               potk += potb;
            }
         } // end if (include)


         shiid = __shfl_sync(ALL_LANES, shiid, ilane + 1);
      }


      if CONSTEXPR (do_g) {
         atomic_add(shgxi[threadIdx.x], gx, shi);
         atomic_add(shgyi[threadIdx.x], gy, shi);
         atomic_add(shgzi[threadIdx.x], gz, shi);
         atomic_add(shtxi[threadIdx.x], &ufld[shi][0]);
         atomic_add(shtyi[threadIdx.x], &ufld[shi][1]);
         atomic_add(shtzi[threadIdx.x], &ufld[shi][2]);
         atomic_add(shdui0[threadIdx.x], &dufld[shi][0]);
         atomic_add(shdui1[threadIdx.x], &dufld[shi][1]);
         atomic_add(shdui2[threadIdx.x], &dufld[shi][2]);
         atomic_add(shdui3[threadIdx.x], &dufld[shi][3]);
         atomic_add(shdui4[threadIdx.x], &dufld[shi][4]);
         atomic_add(shdui5[threadIdx.x], &dufld[shi][5]);
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
      }
      if CONSTEXPR (CFLX) {
         atomic_add(shpoti[threadIdx.x], pot, shi);
         atomic_add(potk, pot, k);
      }
   }
   // */


   //* /
   // block-atoms
   for (int iw = iwarp; iw < niak; iw += nwarp) {
      if CONSTEXPR (do_g) {
         shgxi[threadIdx.x] = 0;
         shgyi[threadIdx.x] = 0;
         shgzi[threadIdx.x] = 0;
         shtxi[threadIdx.x] = 0;
         shtyi[threadIdx.x] = 0;
         shtzi[threadIdx.x] = 0;
         shdui0[threadIdx.x] = 0;
         shdui1[threadIdx.x] = 0;
         shdui2[threadIdx.x] = 0;
         shdui3[threadIdx.x] = 0;
         shdui4[threadIdx.x] = 0;
         shdui5[threadIdx.x] = 0;
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
      }
      if CONSTEXPR (CFLX) {
         shpoti[threadIdx.x] = 0;
         potk = 0;
      }


      int ty = iak[iw];
      int shatomi = ty * WARP_SIZE + ilane;
      int shi = sorted[shatomi].unsorted;
      int atomk = lst[iw * WARP_SIZE + ilane];
      int k = sorted[atomk].unsorted;
      shxi[threadIdx.x] = sorted[shatomi].x;
      shyi[threadIdx.x] = sorted[shatomi].y;
      shzi[threadIdx.x] = sorted[shatomi].z;
      xk = sorted[atomk].x;
      yk = sorted[atomk].y;
      zk = sorted[atomk].z;


      shci[threadIdx.x] = rpole[shi][mpl_pme_0];
      shdix[threadIdx.x] = rpole[shi][mpl_pme_x];
      shdiy[threadIdx.x] = rpole[shi][mpl_pme_y];
      shdiz[threadIdx.x] = rpole[shi][mpl_pme_z];
      shqixx[threadIdx.x] = rpole[shi][mpl_pme_xx];
      shqixy[threadIdx.x] = rpole[shi][mpl_pme_xy];
      shqixz[threadIdx.x] = rpole[shi][mpl_pme_xz];
      shqiyy[threadIdx.x] = rpole[shi][mpl_pme_yy];
      shqiyz[threadIdx.x] = rpole[shi][mpl_pme_yz];
      shqizz[threadIdx.x] = rpole[shi][mpl_pme_zz];
      shuix[threadIdx.x] = uind[shi][0];
      shuiy[threadIdx.x] = uind[shi][1];
      shuiz[threadIdx.x] = uind[shi][2];
      shcorei[threadIdx.x] = pcore[shi];
      shalphai[threadIdx.x] = palpha[shi];
      shvali[threadIdx.x] = pval[shi];
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
         real xi = shxi[klane];
         real yi = shyi[klane];
         real zi = shzi[klane];
         real ci = shci[klane];
         real dix = shdix[klane];
         real diy = shdiy[klane];
         real diz = shdiz[klane];
         real qixx = shqixx[klane];
         real qixy = shqixy[klane];
         real qixz = shqixz[klane];
         real qiyy = shqiyy[klane];
         real qiyz = shqiyz[klane];
         real qizz = shqizz[klane];
         real uix = shuix[klane];
         real uiy = shuiy[klane];
         real uiz = shuiz[klane];
         real corei = shcorei[klane];
         real alphai = shalphai[klane];
         real vali = shvali[klane];


         bool incl = atomk > 0;
         real scaleb = 1;
         real scalec = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;

         real e;
         real pota = 0, potb = 0;
         PairPolarGrad pgrad;
         zero(pgrad);

         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_polar_chgpen<do_e, do_g, ETYP, CFLX>( //
               r2, xr, yr, zr, scaleb, scalec,         //
               ci, dix, diy, diz, corei, vali, alphai, qixx, qixy, qixz, qiyy,
               qiyz, qizz, uix, uiy, uiz, //
               ck, dkx, dky, dkz, corek, valk, alphak, qkxx, qkxy, qkxz, qkyy,
               qkyz, qkzz, ukx, uky, ukz, f, aewald, //
               e, pota, potb, pgrad);

            if CONSTEXPR (do_a)
               if (e != 0)
                  nptl += 1;
            if CONSTEXPR (do_e)
               eptl += cvt_to<ebuf_prec>(e);
            if CONSTEXPR (do_g) {
               shgxi[klane] += pgrad.frcx;
               shgyi[klane] += pgrad.frcy;
               shgzi[klane] += pgrad.frcz;
               gxk -= pgrad.frcx;
               gyk -= pgrad.frcy;
               gzk -= pgrad.frcz;


               shtxi[klane] += pgrad.ufldi[0];
               shtyi[klane] += pgrad.ufldi[1];
               shtzi[klane] += pgrad.ufldi[2];
               txk += pgrad.ufldk[0];
               tyk += pgrad.ufldk[1];
               tzk += pgrad.ufldk[2];


               shdui0[klane] += pgrad.dufldi[0];
               shdui1[klane] += pgrad.dufldi[1];
               shdui2[klane] += pgrad.dufldi[2];
               shdui3[klane] += pgrad.dufldi[3];
               shdui4[klane] += pgrad.dufldi[4];
               shdui5[klane] += pgrad.dufldi[5];
               duk0 += pgrad.dufldk[0];
               duk1 += pgrad.dufldk[1];
               duk2 += pgrad.dufldk[2];
               duk3 += pgrad.dufldk[3];
               duk4 += pgrad.dufldk[4];
               duk5 += pgrad.dufldk[5];
            }
            if CONSTEXPR (do_v) {
               vptlxx += cvt_to<vbuf_prec>(-xr * pgrad.frcx);
               vptlyx += cvt_to<vbuf_prec>(-0.5f *
                                           (yr * pgrad.frcx + xr * pgrad.frcy));
               vptlzx += cvt_to<vbuf_prec>(-0.5f *
                                           (zr * pgrad.frcx + xr * pgrad.frcz));
               vptlyy += cvt_to<vbuf_prec>(-yr * pgrad.frcy);
               vptlzy += cvt_to<vbuf_prec>(-0.5f *
                                           (zr * pgrad.frcy + yr * pgrad.frcz));
               vptlzz += cvt_to<vbuf_prec>(-zr * pgrad.frcz);
            }
            if CONSTEXPR (CFLX) {
               shpoti[klane] += pota;
               potk += potb;
            }
         } // end if (include)
      }


      if CONSTEXPR (do_g) {
         atomic_add(shgxi[threadIdx.x], gx, shi);
         atomic_add(shgyi[threadIdx.x], gy, shi);
         atomic_add(shgzi[threadIdx.x], gz, shi);
         atomic_add(shtxi[threadIdx.x], &ufld[shi][0]);
         atomic_add(shtyi[threadIdx.x], &ufld[shi][1]);
         atomic_add(shtzi[threadIdx.x], &ufld[shi][2]);
         atomic_add(shdui0[threadIdx.x], &dufld[shi][0]);
         atomic_add(shdui1[threadIdx.x], &dufld[shi][1]);
         atomic_add(shdui2[threadIdx.x], &dufld[shi][2]);
         atomic_add(shdui3[threadIdx.x], &dufld[shi][3]);
         atomic_add(shdui4[threadIdx.x], &dufld[shi][4]);
         atomic_add(shdui5[threadIdx.x], &dufld[shi][5]);
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
      }
      if CONSTEXPR (CFLX) {
         atomic_add(shpoti[threadIdx.x], pot, shi);
         atomic_add(potk, pot, k);
      }
   }
   // */


   if CONSTEXPR (do_a) {
      atomic_add(nptl, np, ithread);
   }
   if CONSTEXPR (do_e) {
      atomic_add(eptl, ep, ithread);
   }
   if CONSTEXPR (do_v) {
      atomic_add(vptlxx, vptlyx, vptlzx, vptlyy, vptlzy, vptlzz, vp, ithread);
   }
} // generated by ComplexKernelBuilder (ck.py) 1.5.1


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
   auto bufsize = buffer_size();


   const real f = 0.5f * electric / dielec;
   real aewald = 0;
   if CONSTEXPR (eq<ETYP, EWALD>()) {
      PMEUnit pu = ppme_unit;
      aewald = pu->aewald;
   }


   if CONSTEXPR (do_g)
      darray::zero(g::q0, n, ufld, dufld);

   int ngrid = get_grid_size(BLOCK_DIM);
   epolar_chgpen_cu1<Ver, ETYP, CFLX><<<ngrid, BLOCK_DIM, 0, nonblk>>>(
      st.n, TINKER_IMAGE_ARGS, nep, ep, vir_ep, depx, depy, depz, off,
      st.si1.bit0, nmdwexclude, mdwexclude, mdwexclude_scale, st.x, st.y, st.z,
      st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst, ufld, dufld, uind,
      pot, rpole, pcore, pval, palpha, aewald, f);

   // torque
   if CONSTEXPR (do_g) {
      launch_k1s(nonblk, n, epolar_trq_cu, //
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
