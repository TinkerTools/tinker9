#include "add.h"
#include "cflux.h"
#include "empole_chgpen.h"
#include "empole_chgpen_self.h"
#include "glob.pme.h"
#include "glob.spatial.h"
#include "image.h"
#include "launch.h"
#include "md.h"
#include "pmestuf.h"
#include "seq_bsplgen.h"
#include "seq_pair_mpole_chgpen.h"
#include "seq_triangle.h"
#include "switch.h"
#include "tool/gpu_card.h"


namespace tinker {
template <class Ver, class ETYP, int CFLX>
__global__
void empole_chgpen_cu1(
   int n, TINKER_IMAGE_PARAMS, count_buffer restrict nem,
   energy_buffer restrict em, virial_buffer restrict vem,
   grad_prec* restrict gx, grad_prec* restrict gy, grad_prec* restrict gz,
   real off, const unsigned* restrict minfo, int nexclude,
   const int (*restrict exclude)[2], const real (*restrict exclude_scale)[3],
   const real* restrict x, const real* restrict y, const real* restrict z,
   const Spatial::SortedAtom* restrict sorted, int nakpl,
   const int* restrict iakpl, int niak, const int* restrict iak,
   const int* restrict lst, real* restrict trqx, real* restrict trqy,
   real* restrict trqz, real* restrict pot, const real (*restrict rpole)[10],
   real* restrict pcore, real* restrict pval, const real* restrict palpha,
   real aewald, real f)
{
   constexpr bool do_a = Ver::a;
   constexpr bool do_e = Ver::e;
   constexpr bool do_v = Ver::v;
   constexpr bool do_g = Ver::g;


   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);


   int nemtl;
   if CONSTEXPR (do_a) {
      nemtl = 0;
   }
   using ebuf_prec = energy_buffer_traits::type;
   ebuf_prec emtl;
   if CONSTEXPR (do_e) {
      emtl = 0;
   }
   using vbuf_prec = virial_buffer_traits::type;
   vbuf_prec vemtlxx, vemtlyx, vemtlzx, vemtlyy, vemtlzy, vemtlzz;
   if CONSTEXPR (do_v) {
      vemtlxx = 0;
      vemtlyx = 0;
      vemtlzx = 0;
      vemtlyy = 0;
      vemtlzy = 0;
      vemtlzz = 0;
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
   __shared__ real shpoti[BLOCK_DIM];
   real gxk;
   real gyk;
   real gzk;
   real txk;
   real tyk;
   real tzk;
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
         gxk = 0;
         gyk = 0;
         gzk = 0;
         txk = 0;
         tyk = 0;
         tzk = 0;
      }
      if CONSTEXPR (CFLX) {
         shpoti[threadIdx.x] = 0;
         potk = 0;
      }



      int shi = exclude[ii][0];
      int k = exclude[ii][1];
      real scalea = exclude_scale[ii][0];


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
      corek = pcore[k];
      alphak = palpha[k];
      valk = pval[k];


      int klane = threadIdx.x;
      constexpr bool incl = true;
      real xr = xk - xi;
      real yr = yk - yi;
      real zr = zk - zi;

      real e;
      real pota, potb;
      PairMPoleGrad pgrad;
      zero(pgrad);

      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         pair_mpole_chgpen<do_e, do_g, ETYP, CFLX>(
            r2, xr, yr, zr, scalea,                 //
            ci, dix, diy, diz, corei, vali, alphai, //
            qixx, qixy, qixz, qiyy, qiyz, qizz,     //
            ck, dkx, dky, dkz, corek, valk, alphak, //
            qkxx, qkxy, qkxz, qkyy, qkyz, qkzz,     //
            f, aewald, e, pota, potb, pgrad);
         
         
         if CONSTEXPR (do_a and scalea != 0)
            if (e != 0)
               nemtl += 1;
         if CONSTEXPR (do_e)
            emtl += cvt_to<ebuf_prec>(e);
         if CONSTEXPR (do_g) {
            shgxi[klane] += pgrad.frcx;
            shgyi[klane] += pgrad.frcy;
            shgzi[klane] += pgrad.frcz;
            gxk -= pgrad.frcx;
            gyk -= pgrad.frcy;
            gzk -= pgrad.frcz;

            shtxi[klane] += pgrad.ttmi[0];
            shtyi[klane] += pgrad.ttmi[1];
            shtzi[klane] += pgrad.ttmi[2];
            txk += pgrad.ttmk[0];
            tyk += pgrad.ttmk[1];
            tzk += pgrad.ttmk[2];
         }
         if CONSTEXPR (do_v) {
            vemtlxx += cvt_to<vbuf_prec>(-xr * pgrad.frcx);
            vemtlyx +=
               cvt_to<vbuf_prec>(-0.5f * (yr * pgrad.frcx + xr * pgrad.frcy));
            vemtlzx +=
               cvt_to<vbuf_prec>(-0.5f * (zr * pgrad.frcx + xr * pgrad.frcz));
            vemtlyy += cvt_to<vbuf_prec>(-yr * pgrad.frcy);
            vemtlzy +=
               cvt_to<vbuf_prec>(-0.5f * (zr * pgrad.frcy + yr * pgrad.frcz));
            vemtlzz += cvt_to<vbuf_prec>(-zr * pgrad.frcz);
         }
         if CONSTEXPR (CFLX) {
            shpoti[klane] += pota;
            potk += potb;
            printf("%2d %2d %15.8e %15.8e\n", klane+1,k+1,pota,potb);
         }
      } // end if (include)


      if CONSTEXPR (do_g) {
         atomic_add(shgxi[threadIdx.x], gx, shi);
         atomic_add(shgyi[threadIdx.x], gy, shi);
         atomic_add(shgzi[threadIdx.x], gz, shi);
         atomic_add(shtxi[threadIdx.x], trqx, shi);
         atomic_add(shtyi[threadIdx.x], trqy, shi);
         atomic_add(shtzi[threadIdx.x], trqz, shi);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
         atomic_add(txk, trqx, k);
         atomic_add(tyk, trqy, k);
         atomic_add(tzk, trqz, k);
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
         gxk = 0;
         gyk = 0;
         gzk = 0;
         txk = 0;
         tyk = 0;
         tzk = 0;
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
      corek = pcore[k];
      alphak = palpha[k];
      valk = pval[k];


      unsigned int minfo0 = minfo[iw * WARP_SIZE + ilane];


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
         real corei = shcorei[klane];
         real alphai = shalphai[klane];
         real vali = shvali[klane];


         bool incl = iid < kid and kid < n;
         incl = incl and (minfo0 & srcmask) == 0;
         real scalea = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;

         real e;
         real pota, potb;
         PairMPoleGrad pgrad;
         zero(pgrad);

         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_mpole_chgpen<do_e, do_g, ETYP, CFLX>(
               r2, xr, yr, zr, scalea,                 //
               ci, dix, diy, diz, corei, vali, alphai, //
               qixx, qixy, qixz, qiyy, qiyz, qizz,     //
               ck, dkx, dky, dkz, corek, valk, alphak, //
               qkxx, qkxy, qkxz, qkyy, qkyz, qkzz,     //
               f, aewald, e, pota, potb, pgrad);


            if CONSTEXPR (do_a)
               if (e != 0)
                  nemtl += 1;
            if CONSTEXPR (do_e)
               emtl += cvt_to<ebuf_prec>(e);
            if CONSTEXPR (do_g) {
               shgxi[klane] += pgrad.frcx;
               shgyi[klane] += pgrad.frcy;
               shgzi[klane] += pgrad.frcz;
               gxk -= pgrad.frcx;
               gyk -= pgrad.frcy;
               gzk -= pgrad.frcz;

               shtxi[klane] += pgrad.ttmi[0];
               shtyi[klane] += pgrad.ttmi[1];
               shtzi[klane] += pgrad.ttmi[2];
               txk += pgrad.ttmk[0];
               tyk += pgrad.ttmk[1];
               tzk += pgrad.ttmk[2];
            }
            if CONSTEXPR (do_v) {
               vemtlxx += cvt_to<vbuf_prec>(-xr * pgrad.frcx);
               vemtlyx += cvt_to<vbuf_prec>(
                  -0.5f * (yr * pgrad.frcx + xr * pgrad.frcy));
               vemtlzx += cvt_to<vbuf_prec>(
                  -0.5f * (zr * pgrad.frcx + xr * pgrad.frcz));
               vemtlyy += cvt_to<vbuf_prec>(-yr * pgrad.frcy);
               vemtlzy += cvt_to<vbuf_prec>(
                  -0.5f * (zr * pgrad.frcy + yr * pgrad.frcz));
               vemtlzz += cvt_to<vbuf_prec>(-zr * pgrad.frcz);
            }
            if CONSTEXPR (CFLX) {
               shpoti[klane] += pota;
               potk += potb;
               printf("%2d %2d %15.8e %15.8e\n", klane+1,k+1,pota,potb);
            }
         } // end if (include)


         shiid = __shfl_sync(ALL_LANES, shiid, ilane + 1);
      }


      if CONSTEXPR (do_g) {
         atomic_add(shgxi[threadIdx.x], gx, shi);
         atomic_add(shgyi[threadIdx.x], gy, shi);
         atomic_add(shgzi[threadIdx.x], gz, shi);
         atomic_add(shtxi[threadIdx.x], trqx, shi);
         atomic_add(shtyi[threadIdx.x], trqy, shi);
         atomic_add(shtzi[threadIdx.x], trqz, shi);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
         atomic_add(txk, trqx, k);
         atomic_add(tyk, trqy, k);
         atomic_add(tzk, trqz, k);
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
         gxk = 0;
         gyk = 0;
         gzk = 0;
         txk = 0;
         tyk = 0;
         tzk = 0;
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
         real corei = shcorei[klane];
         real alphai = shalphai[klane];
         real vali = shvali[klane];


         bool incl = atomk > 0;
         real scalea = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;

         real e;
         real pota, potb;
         PairMPoleGrad pgrad;
         zero(pgrad);

         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_mpole_chgpen<do_e, do_g, ETYP, CFLX>(
               r2, xr, yr, zr, scalea,                 //
               ci, dix, diy, diz, corei, vali, alphai, //
               qixx, qixy, qixz, qiyy, qiyz, qizz,     //
               ck, dkx, dky, dkz, corek, valk, alphak, //
               qkxx, qkxy, qkxz, qkyy, qkyz, qkzz,     //
               f, aewald, e, pota, potb, pgrad);


            if CONSTEXPR (do_a)
               if (e != 0)
                  nemtl += 1;
            if CONSTEXPR (do_e)
               emtl += cvt_to<ebuf_prec>(e);
            if CONSTEXPR (do_g) {
               shgxi[klane] += pgrad.frcx;
               shgyi[klane] += pgrad.frcy;
               shgzi[klane] += pgrad.frcz;
               gxk -= pgrad.frcx;
               gyk -= pgrad.frcy;
               gzk -= pgrad.frcz;

               shtxi[klane] += pgrad.ttmi[0];
               shtyi[klane] += pgrad.ttmi[1];
               shtzi[klane] += pgrad.ttmi[2];
               txk += pgrad.ttmk[0];
               tyk += pgrad.ttmk[1];
               tzk += pgrad.ttmk[2];
            }
            if CONSTEXPR (do_v) {
               vemtlxx += cvt_to<vbuf_prec>(-xr * pgrad.frcx);
               vemtlyx += cvt_to<vbuf_prec>(
                  -0.5f * (yr * pgrad.frcx + xr * pgrad.frcy));
               vemtlzx += cvt_to<vbuf_prec>(
                  -0.5f * (zr * pgrad.frcx + xr * pgrad.frcz));
               vemtlyy += cvt_to<vbuf_prec>(-yr * pgrad.frcy);
               vemtlzy += cvt_to<vbuf_prec>(
                  -0.5f * (zr * pgrad.frcy + yr * pgrad.frcz));
               vemtlzz += cvt_to<vbuf_prec>(-zr * pgrad.frcz);
            }
            if CONSTEXPR (CFLX) {
               shpoti[klane] += pota;
               potk += potb;
               printf("%2d %2d %15.8e %15.8e\n", klane+1,k+1,pota,potb);
            }
         } // end if (include)
      }


      if CONSTEXPR (do_g) {
         atomic_add(shgxi[threadIdx.x], gx, shi);
         atomic_add(shgyi[threadIdx.x], gy, shi);
         atomic_add(shgzi[threadIdx.x], gz, shi);
         atomic_add(shtxi[threadIdx.x], trqx, shi);
         atomic_add(shtyi[threadIdx.x], trqy, shi);
         atomic_add(shtzi[threadIdx.x], trqz, shi);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
         atomic_add(txk, trqx, k);
         atomic_add(tyk, trqy, k);
         atomic_add(tzk, trqz, k);
      }
      if CONSTEXPR (CFLX) {
        atomic_add(shpoti[threadIdx.x], pot, shi);
        atomic_add(potk, pot, k);
      }
   }
   // */


   if CONSTEXPR (do_a) {
      atomic_add(nemtl, nem, ithread);
   }
   if CONSTEXPR (do_e) {
      atomic_add(emtl, em, ithread);
   }
   if CONSTEXPR (do_v) {
      atomic_add(vemtlxx, vemtlyx, vemtlzx, vemtlyy, vemtlzy, vemtlzz, vem,
                 ithread);
   }
} // generated by ComplexKernelBuilder (ck.py) 1.5.1


template <class Ver, class ETYP, int CFLX>
void empole_chgpen_cu()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;


   // const auto& st = *mspatial_unit;
   const auto& st = *mspatial_v2_unit;
   real off;
   if CONSTEXPR (eq<ETYP, EWALD>())
      off = switch_off(switch_ewald);
   else
      off = switch_off(switch_mpole);
   // const real off2 = off * off;
   auto bufsize = buffer_size();


   const real f = electric / dielec;
   real aewald = 0;
   if CONSTEXPR (eq<ETYP, EWALD>()) {
      PMEUnit pu = epme_unit;
      aewald = pu->aewald;


      
      launch_k1s(nonblk, n, empole_chgpen_self_cu<do_a, do_e, CFLX>, //
                  bufsize, nem, em, rpole, pot, n, f, aewald);
   }


   int ngrid = get_grid_size(BLOCK_DIM);
   empole_chgpen_cu1<Ver, ETYP, CFLX><<<ngrid, BLOCK_DIM, 0, nonblk>>>(
      st.n, TINKER_IMAGE_ARGS, nem, em, vir_em, demx, demy, demz, off,
      st.si1.bit0, nmdwexclude, mdwexclude, mdwexclude_scale, st.x, st.y, st.z,
      st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst, trqx, trqy, trqz,
      pot, rpole, pcore, pval, palpha, aewald, f);
}


void empole_chgpen_nonewald_cu(int vers, int use_cf)
{
   if (use_cf) {
      if (vers == calc::v0) {
         empole_chgpen_cu<calc::V0, NON_EWALD, 1>();
      } else if (vers == calc::v1) {
         empole_chgpen_cu<calc::V1, NON_EWALD, 1>();
      } else if (vers == calc::v3) {
         empole_chgpen_cu<calc::V3, NON_EWALD, 1>();
      } else if (vers == calc::v4) {
         empole_chgpen_cu<calc::V4, NON_EWALD, 1>();
      } else if (vers == calc::v5) {
         empole_chgpen_cu<calc::V5, NON_EWALD, 1>();
      } else if (vers == calc::v6) {
         empole_chgpen_cu<calc::V6, NON_EWALD, 1>();
      }
   } else {
      if (vers == calc::v0) {
         empole_chgpen_cu<calc::V0, NON_EWALD, 0>();
      } else if (vers == calc::v1) {
         empole_chgpen_cu<calc::V1, NON_EWALD, 0>();
      } else if (vers == calc::v3) {
         empole_chgpen_cu<calc::V3, NON_EWALD, 0>();
      } else if (vers == calc::v4) {
         empole_chgpen_cu<calc::V4, NON_EWALD, 0>();
      } else if (vers == calc::v5) {
         empole_chgpen_cu<calc::V5, NON_EWALD, 0>();
      } else if (vers == calc::v6) {
         empole_chgpen_cu<calc::V6, NON_EWALD, 0>();
      }
   }
}


void empole_chgpen_ewald_real_self_cu(int vers, int use_cf)
{
   if (use_cf) {
      if (vers == calc::v0) {
         empole_chgpen_cu<calc::V0, EWALD, 1>();
      } else if (vers == calc::v1) {
         empole_chgpen_cu<calc::V1, EWALD, 1>();
      } else if (vers == calc::v3) {
         empole_chgpen_cu<calc::V3, EWALD, 1>();
      } else if (vers == calc::v4) {
         empole_chgpen_cu<calc::V4, EWALD, 1>();
      } else if (vers == calc::v5) {
         empole_chgpen_cu<calc::V5, EWALD, 1>();
      } else if (vers == calc::v6) {
         empole_chgpen_cu<calc::V6, EWALD, 1>();
      }
   } else {
      if (vers == calc::v0) {
         empole_chgpen_cu<calc::V0, EWALD, 0>();
      } else if (vers == calc::v1) {
         empole_chgpen_cu<calc::V1, EWALD, 0>();
      } else if (vers == calc::v3) {
         empole_chgpen_cu<calc::V3, EWALD, 0>();
      } else if (vers == calc::v4) {
         empole_chgpen_cu<calc::V4, EWALD, 0>();
      } else if (vers == calc::v5) {
         empole_chgpen_cu<calc::V5, EWALD, 0>();
      } else if (vers == calc::v6) {
         empole_chgpen_cu<calc::V6, EWALD, 0>();
      }
   }
}
}
