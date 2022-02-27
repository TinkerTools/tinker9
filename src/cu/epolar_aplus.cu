#include "add.h"
#include "cflux.h"
#include "epolar_aplus.h"
#include "epolar_trq.h"
#include "glob.mplar.h"
#include "glob.spatial.h"
#include "image.h"
#include "launch.h"
#include "md.h"
#include "pme.h"
#include "seq_pair_polar_aplus.h"
#include "seq_triangle.h"
#include "switch.h"
#include "tool/gpu_card.h"


namespace tinker {
// ck.py Version 2.0.3
template <class Ver, class ETYP, bool CFLX>
__global__
void epolar_aplus_cu1(int n, TINKER_IMAGE_PARAMS, count_buffer restrict nep,
                energy_buffer restrict ep, virial_buffer restrict vep,
                grad_prec* restrict gx, grad_prec* restrict gy,
                grad_prec* restrict gz, real off,
                const unsigned* restrict mdpuinfo, int nexclude,
                const int (*restrict exclude)[2],
                const real (*restrict exclude_scale)[4], const real* restrict x,
                const real* restrict y, const real* restrict z,
                const Spatial::SortedAtom* restrict sorted, int nakpl,
                const int* restrict iakpl, int niak, const int* restrict iak,
                const int* restrict lst, real (*restrict ufld)[3],
                real (*restrict dufld)[6], const real (*restrict uind)[3], real* restrict pot,
                const real (*restrict rpole)[10], const real* restrict pdamp, const real* restrict thole, 
	              const real* restrict dirdamp, real aewald, real f)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);


   int neptl;
   if CONSTEXPR (do_a) {
      neptl = 0;
   }
   using ebuf_prec = energy_buffer_traits::type;
   ebuf_prec eptl;
   if CONSTEXPR (do_e) {
      eptl = 0;
   }
   using vbuf_prec = virial_buffer_traits::type;
   vbuf_prec veptlxx, veptlyx, veptlzx, veptlyy, veptlzy, veptlzz;
   if CONSTEXPR (do_v) {
      veptlxx = 0;
      veptlyx = 0;
      veptlzx = 0;
      veptlyy = 0;
      veptlzy = 0;
      veptlzz = 0;
   }
   __shared__ real xi[BLOCK_DIM];
   __shared__ real yi[BLOCK_DIM];
   __shared__ real zi[BLOCK_DIM];
   real xk;
   real yk;
   real zk;
   __shared__ real frcxi[BLOCK_DIM];
   __shared__ real frcyi[BLOCK_DIM];
   __shared__ real frczi[BLOCK_DIM];
   __shared__ real ufld0i[BLOCK_DIM];
   __shared__ real ufld1i[BLOCK_DIM];
   __shared__ real ufld2i[BLOCK_DIM];
   __shared__ real dufld0i[BLOCK_DIM];
   __shared__ real dufld1i[BLOCK_DIM];
   __shared__ real dufld2i[BLOCK_DIM];
   __shared__ real dufld3i[BLOCK_DIM];
   __shared__ real dufld4i[BLOCK_DIM];
   __shared__ real dufld5i[BLOCK_DIM];
   __shared__ real poti[BLOCK_DIM];
   real frcxk;
   real frcyk;
   real frczk;
   real ufld0k;
   real ufld1k;
   real ufld2k;
   real dufld0k;
   real dufld1k;
   real dufld2k;
   real dufld3k;
   real dufld4k;
   real dufld5k;
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
   __shared__ real pdi[BLOCK_DIM];
   __shared__ real pti[BLOCK_DIM];
   __shared__ real ddi[BLOCK_DIM];
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
   real pdk;
   real ptk;
   real ddk;


   //* /
   for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
      const int klane = threadIdx.x;
      if CONSTEXPR (do_g) {
         frcxi[threadIdx.x] = 0;
         frcyi[threadIdx.x] = 0;
         frczi[threadIdx.x] = 0;
         ufld0i[threadIdx.x] = 0;
         ufld1i[threadIdx.x] = 0;
         ufld2i[threadIdx.x] = 0;
         dufld0i[threadIdx.x] = 0;
         dufld1i[threadIdx.x] = 0;
         dufld2i[threadIdx.x] = 0;
         dufld3i[threadIdx.x] = 0;
         dufld4i[threadIdx.x] = 0;
         dufld5i[threadIdx.x] = 0;
         if CONSTEXPR (CFLX)
            poti[threadIdx.x] = 0;
         frcxk = 0;
         frcyk = 0;
         frczk = 0;
         ufld0k = 0;
         ufld1k = 0;
         ufld2k = 0;
         dufld0k = 0;
         dufld1k = 0;
         dufld2k = 0;
         dufld3k = 0;
         dufld4k = 0;
         dufld5k = 0;
         if CONSTEXPR (CFLX)
            potk = 0;
      }


      int i = exclude[ii][0];
      int k = exclude[ii][1];
      real scaleb = exclude_scale[ii][1]; //d
      real scaled = exclude_scale[ii][3]; //u

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
      pdi[klane] = pdamp[i];
      pti[klane] = thole[i];
      ddi[klane] = dirdamp[i];
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
      pdk = pdamp[k];
      ptk = thole[k];
      ddk = dirdamp[k];


      constexpr bool incl = true;
      real xr = xk - xi[klane];
      real yr = yk - yi[klane];
      real zr = zk - zi[klane];
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         real e, vxx, vyx, vzx, vyy, vzy, vzz;
         real e1, vxx1, vyx1, vzx1, vyy1, vzy1, vzz1;
      	 real pota, potb;
      	 real pota1, potb1;
         pair_polar_aplus_v2<Ver, ETYP, CFLX>( //
            r2, xr, yr, zr, 1, 1, //
            ci[klane], dix[klane], diy[klane], diz[klane], qixx[klane], //
            qixy[klane], qixz[klane], qiyy[klane], qiyz[klane], qizz[klane], //
            uix[klane], uiy[klane], uiz[klane], //
            pdi[klane], pti[klane], ddi[klane], //
            ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukx,  //
            uky, ukz, pdk, ptk, ddk,  //
            f, aewald, //
            frcxi[klane], frcyi[klane], frczi[klane], frcxk, frcyk, frczk,
            ufld0i[klane], ufld1i[klane], ufld2i[klane], ufld0k, ufld1k,
            ufld2k, //
            dufld0i[klane], dufld1i[klane], dufld2i[klane], dufld3i[klane],
            dufld4i[klane], dufld5i[klane], dufld0k, dufld1k, dufld2k, dufld3k,
            dufld4k, dufld5k, //
            e1, vxx1, vyx1, vzx1, vyy1, vzy1, vzz1, pota1, potb1);
         pair_polar_aplus_v2<Ver, NON_EWALD, CFLX>( //
            r2, xr, yr, zr, scaleb - 1, scaled - 1, //
            ci[klane], dix[klane], diy[klane], diz[klane], qixx[klane], //
            qixy[klane], qixz[klane], qiyy[klane], qiyz[klane], qizz[klane], //
            uix[klane], uiy[klane], uiz[klane], //
            pdi[klane], pti[klane], ddi[klane], //
            ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukx,  //
            uky, ukz, pdk, ptk, ddk,  //
            f, aewald, //
            frcxi[klane], frcyi[klane], frczi[klane], frcxk, frcyk, frczk,
            ufld0i[klane], ufld1i[klane], ufld2i[klane], ufld0k, ufld1k,
            ufld2k, //
            dufld0i[klane], dufld1i[klane], dufld2i[klane], dufld3i[klane],
            dufld4i[klane], dufld5i[klane], dufld0k, dufld1k, dufld2k, dufld3k,
            dufld4k, dufld5k, //
            e, vxx, vyx, vzx, vyy, vzy, vzz, pota, potb);
         if CONSTEXPR (do_e) {
            e = e + e1;
            eptl += cvt_to<ebuf_prec>(e);
            if CONSTEXPR (do_a) {
               if (e != 0 and scaleb != 0)
               neptl += 1;
            }
         }
         if CONSTEXPR (do_v) {
            veptlxx += cvt_to<vbuf_prec>(vxx + vxx1);
            veptlyx += cvt_to<vbuf_prec>(vyx + vyx1);
            veptlzx += cvt_to<vbuf_prec>(vzx + vzx1);
            veptlyy += cvt_to<vbuf_prec>(vyy + vyy1);
            veptlzy += cvt_to<vbuf_prec>(vzy + vzy1);
            veptlzz += cvt_to<vbuf_prec>(vzz + vzz1);
         }
         if CONSTEXPR (CFLX) {
            poti[klane] += (pota + pota1);
            potk += (potb + potb1);
         }
      } // end if (include)


      if CONSTEXPR (do_g) {
         atomic_add(frcxi[threadIdx.x], gx, i);
         atomic_add(frcyi[threadIdx.x], gy, i);
         atomic_add(frczi[threadIdx.x], gz, i);
         atomic_add(ufld0i[threadIdx.x], &ufld[i][0]);
         atomic_add(ufld1i[threadIdx.x], &ufld[i][1]);
         atomic_add(ufld2i[threadIdx.x], &ufld[i][2]);
         atomic_add(dufld0i[threadIdx.x], &dufld[i][0]);
         atomic_add(dufld1i[threadIdx.x], &dufld[i][1]);
         atomic_add(dufld2i[threadIdx.x], &dufld[i][2]);
         atomic_add(dufld3i[threadIdx.x], &dufld[i][3]);
         atomic_add(dufld4i[threadIdx.x], &dufld[i][4]);
         atomic_add(dufld5i[threadIdx.x], &dufld[i][5]);
         if CONSTEXPR (CFLX)
            atomic_add(poti[threadIdx.x], pot, i);
         atomic_add(frcxk, gx, k);
         atomic_add(frcyk, gy, k);
         atomic_add(frczk, gz, k);
         atomic_add(ufld0k, &ufld[k][0]);
         atomic_add(ufld1k, &ufld[k][1]);
         atomic_add(ufld2k, &ufld[k][2]);
         atomic_add(dufld0k, &dufld[k][0]);
         atomic_add(dufld1k, &dufld[k][1]);
         atomic_add(dufld2k, &dufld[k][2]);
         atomic_add(dufld3k, &dufld[k][3]);
         atomic_add(dufld4k, &dufld[k][4]);
         atomic_add(dufld5k, &dufld[k][5]);
         if CONSTEXPR (CFLX)
            atomic_add(potk, pot, k);
      }
   }
   // */


   for (int iw = iwarp; iw < nakpl; iw += nwarp) {
      if CONSTEXPR (do_g) {
         frcxi[threadIdx.x] = 0;
         frcyi[threadIdx.x] = 0;
         frczi[threadIdx.x] = 0;
         ufld0i[threadIdx.x] = 0;
         ufld1i[threadIdx.x] = 0;
         ufld2i[threadIdx.x] = 0;
         dufld0i[threadIdx.x] = 0;
         dufld1i[threadIdx.x] = 0;
         dufld2i[threadIdx.x] = 0;
         dufld3i[threadIdx.x] = 0;
         dufld4i[threadIdx.x] = 0;
         dufld5i[threadIdx.x] = 0;
         if CONSTEXPR (CFLX)
            poti[threadIdx.x] = 0;
         frcxk = 0;
         frcyk = 0;
         frczk = 0;
         ufld0k = 0;
         ufld1k = 0;
         ufld2k = 0;
         dufld0k = 0;
         dufld1k = 0;
         dufld2k = 0;
         dufld3k = 0;
         dufld4k = 0;
         dufld5k = 0;
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
      pdi[threadIdx.x] = pdamp[i];
      pti[threadIdx.x] = thole[i];
      ddi[threadIdx.x] = dirdamp[i];
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
      pdk = pdamp[k];
      ptk = thole[k];
      ddk = dirdamp[k];


      unsigned int mdpuinfo0 = mdpuinfo[iw * WARP_SIZE + ilane];
      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         bool incl = iid < kid and kid < n;
         int srcmask = 1 << srclane;
         incl = incl and (mdpuinfo0 & srcmask) == 0;
         real xr = xk - xi[klane];
         real yr = yk - yi[klane];
         real zr = zk - zi[klane];
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real e, vxx, vyx, vzx, vyy, vzy, vzz;
      	 		real pota, potb;
            pair_polar_aplus_v2<Ver, NON_EWALD, CFLX>( //
               r2, xr, yr, zr, 1, 1, 
               ci[klane], dix[klane], diy[klane], diz[klane], qixx[klane], //
               qixy[klane], qixz[klane], qiyy[klane], qiyz[klane], qizz[klane], //
               uix[klane], uiy[klane], uiz[klane], //
               pdi[klane], pti[klane], ddi[klane], //
               ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukx,  //
               uky, ukz, pdk, ptk, ddk,  //
               f, aewald, //
               frcxi[klane], frcyi[klane], frczi[klane], frcxk, frcyk, frczk,
               ufld0i[klane], ufld1i[klane], ufld2i[klane], ufld0k, ufld1k,
               ufld2k, //
               dufld0i[klane], dufld1i[klane], dufld2i[klane], dufld3i[klane],
               dufld4i[klane], dufld5i[klane], dufld0k, dufld1k, dufld2k,
               dufld3k, dufld4k, dufld5k, //
               e, vxx, vyx, vzx, vyy, vzy, vzz, pota, potb);
            if CONSTEXPR (do_e) {
               eptl += cvt_to<ebuf_prec>(e);
               if CONSTEXPR (do_a) {
                  if (e != 0)
                     neptl += 1;
               }
            }
            if CONSTEXPR (do_v) {
               veptlxx += cvt_to<vbuf_prec>(vxx);
               veptlyx += cvt_to<vbuf_prec>(vyx);
               veptlzx += cvt_to<vbuf_prec>(vzx);
               veptlyy += cvt_to<vbuf_prec>(vyy);
               veptlzy += cvt_to<vbuf_prec>(vzy);
               veptlzz += cvt_to<vbuf_prec>(vzz);
            }
            if CONSTEXPR (CFLX) {
               poti[klane] += pota;
               potk += potb;
            }
         } // end if (include)


         iid = __shfl_sync(ALL_LANES, iid, ilane + 1);
      }


      if CONSTEXPR (do_g) {
         atomic_add(frcxi[threadIdx.x], gx, i);
         atomic_add(frcyi[threadIdx.x], gy, i);
         atomic_add(frczi[threadIdx.x], gz, i);
         atomic_add(ufld0i[threadIdx.x], &ufld[i][0]);
         atomic_add(ufld1i[threadIdx.x], &ufld[i][1]);
         atomic_add(ufld2i[threadIdx.x], &ufld[i][2]);
         atomic_add(dufld0i[threadIdx.x], &dufld[i][0]);
         atomic_add(dufld1i[threadIdx.x], &dufld[i][1]);
         atomic_add(dufld2i[threadIdx.x], &dufld[i][2]);
         atomic_add(dufld3i[threadIdx.x], &dufld[i][3]);
         atomic_add(dufld4i[threadIdx.x], &dufld[i][4]);
         atomic_add(dufld5i[threadIdx.x], &dufld[i][5]);
         if CONSTEXPR (CFLX)
            atomic_add(poti[threadIdx.x], pot, i);
         atomic_add(frcxk, gx, k);
         atomic_add(frcyk, gy, k);
         atomic_add(frczk, gz, k);
         atomic_add(ufld0k, &ufld[k][0]);
         atomic_add(ufld1k, &ufld[k][1]);
         atomic_add(ufld2k, &ufld[k][2]);
         atomic_add(dufld0k, &dufld[k][0]);
         atomic_add(dufld1k, &dufld[k][1]);
         atomic_add(dufld2k, &dufld[k][2]);
         atomic_add(dufld3k, &dufld[k][3]);
         atomic_add(dufld4k, &dufld[k][4]);
         atomic_add(dufld5k, &dufld[k][5]);
         if CONSTEXPR (CFLX)
            atomic_add(potk, pot, k);
      }
   }


   for (int iw = iwarp; iw < niak; iw += nwarp) {
      if CONSTEXPR (do_g) {
         frcxi[threadIdx.x] = 0;
         frcyi[threadIdx.x] = 0;
         frczi[threadIdx.x] = 0;
         ufld0i[threadIdx.x] = 0;
         ufld1i[threadIdx.x] = 0;
         ufld2i[threadIdx.x] = 0;
         dufld0i[threadIdx.x] = 0;
         dufld1i[threadIdx.x] = 0;
         dufld2i[threadIdx.x] = 0;
         dufld3i[threadIdx.x] = 0;
         dufld4i[threadIdx.x] = 0;
         dufld5i[threadIdx.x] = 0;
         frcxk = 0;
         frcyk = 0;
         frczk = 0;
         ufld0k = 0;
         if CONSTEXPR (CFLX)
            poti[threadIdx.x] = 0;
         frcxk = 0;
         frcyk = 0;
         frczk = 0;
         ufld0k = 0;
         ufld1k = 0;
         ufld2k = 0;
         dufld0k = 0;
         dufld1k = 0;
         dufld2k = 0;
         dufld3k = 0;
         dufld4k = 0;
         dufld5k = 0;
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
      pdi[threadIdx.x] = pdamp[i];
      pti[threadIdx.x] = thole[i];
      ddi[threadIdx.x] = dirdamp[i];
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
      pdk = pdamp[k];
      ptk = thole[k];
      ddk = dirdamp[k];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         bool incl = atomk > 0;
         real xr = xk - xi[klane];
         real yr = yk - yi[klane];
         real zr = zk - zi[klane];
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real e, vxx, vyx, vzx, vyy, vzy, vzz;
      	 		real pota, potb;
            pair_polar_aplus_v2<Ver, NON_EWALD, CFLX>( //
               r2, xr, yr, zr, 1, 1, 
               ci[klane], dix[klane], diy[klane], diz[klane], qixx[klane], //
               qixy[klane], qixz[klane], qiyy[klane], qiyz[klane], qizz[klane], //
               uix[klane], uiy[klane], uiz[klane], //
               pdi[klane], pti[klane], ddi[klane], //
               ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukx,  //
               uky, ukz, pdk, ptk, ddk,  //
               f, aewald, //
               frcxi[klane], frcyi[klane], frczi[klane], frcxk, frcyk, frczk,
               ufld0i[klane], ufld1i[klane], ufld2i[klane], ufld0k, ufld1k,
               ufld2k, //
               dufld0i[klane], dufld1i[klane], dufld2i[klane], dufld3i[klane],
               dufld4i[klane], dufld5i[klane], dufld0k, dufld1k, dufld2k,
               dufld3k, dufld4k, dufld5k, //
               e, vxx, vyx, vzx, vyy, vzy, vzz, pota, potb);
						if CONSTEXPR (do_e) {
               eptl += cvt_to<ebuf_prec>(e);
               if CONSTEXPR (do_a) {
                  if (e != 0)
                     neptl += 1;
               }
            }
            if CONSTEXPR (do_v) {
               veptlxx += cvt_to<vbuf_prec>(vxx);
               veptlyx += cvt_to<vbuf_prec>(vyx);
               veptlzx += cvt_to<vbuf_prec>(vzx);
               veptlyy += cvt_to<vbuf_prec>(vyy);
               veptlzy += cvt_to<vbuf_prec>(vzy);
               veptlzz += cvt_to<vbuf_prec>(vzz);
            }
            if CONSTEXPR (CFLX) {
               poti[klane] += pota;
               potk += potb;
            }
         } // end if (include)
      }


      if CONSTEXPR (do_g) {
         atomic_add(frcxi[threadIdx.x], gx, i);
         atomic_add(frcyi[threadIdx.x], gy, i);
         atomic_add(frczi[threadIdx.x], gz, i);
         atomic_add(ufld0i[threadIdx.x], &ufld[i][0]);
         atomic_add(ufld1i[threadIdx.x], &ufld[i][1]);
         atomic_add(ufld2i[threadIdx.x], &ufld[i][2]);
         atomic_add(dufld0i[threadIdx.x], &dufld[i][0]);
         atomic_add(dufld1i[threadIdx.x], &dufld[i][1]);
         atomic_add(dufld2i[threadIdx.x], &dufld[i][2]);
         atomic_add(dufld3i[threadIdx.x], &dufld[i][3]);
         atomic_add(dufld4i[threadIdx.x], &dufld[i][4]);
         atomic_add(dufld5i[threadIdx.x], &dufld[i][5]);
         if CONSTEXPR (CFLX)
            atomic_add(poti[threadIdx.x], pot, i);
         atomic_add(frcxk, gx, k);
         atomic_add(frcyk, gy, k);
         atomic_add(frczk, gz, k);
         atomic_add(ufld0k, &ufld[k][0]);
         atomic_add(ufld1k, &ufld[k][1]);
         atomic_add(ufld2k, &ufld[k][2]);
         atomic_add(dufld0k, &dufld[k][0]);
         atomic_add(dufld1k, &dufld[k][1]);
         atomic_add(dufld2k, &dufld[k][2]);
         atomic_add(dufld3k, &dufld[k][3]);
         atomic_add(dufld4k, &dufld[k][4]);
         atomic_add(dufld5k, &dufld[k][5]);
         if CONSTEXPR (CFLX)
            atomic_add(potk, pot, k);
      }
   }


   if CONSTEXPR (do_a) {
      atomic_add(neptl, nep, ithread);
   }
   if CONSTEXPR (do_e) {
      atomic_add(eptl, ep, ithread);
   }
   if CONSTEXPR (do_v) {
      atomic_add(veptlxx, veptlyx, veptlzx, veptlyy, veptlzy, veptlzz, vep,
                 ithread);
   }
}


template <class Ver, class ETYP, int CFLX>
void epolar_aplus_cu(const real (*uind)[3])
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


   if CONSTEXPR (do_g) {
      darray::zero(g::q0, n, ufld, dufld);
   }
   int ngrid = get_grid_size(BLOCK_DIM);
   epolar_aplus_cu1<Ver, ETYP, CFLX><<<ngrid, BLOCK_DIM, 0, g::s0>>>(
      st.n, TINKER_IMAGE_ARGS, nep, ep, vir_ep, depx, depy, depz, off,
      st.si1.bit0, nmdpuexclude, mdpuexclude, mdpuexclude_scale, st.x, st.y, st.z,
      st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst, ufld, dufld, uind,
      pot, rpole, pdamp, thole, dirdamp, aewald, f);


   // torque
   if CONSTEXPR (do_g) {
      launch_k1s(g::s0, n, epolar_trq_cu, //
                 trqx, trqy, trqz, n, rpole, ufld, dufld);
   }
}


void epolar_aplus_nonewald_cu(int vers, int use_cf, const real (*uind)[3])
{
   if (use_cf) {
      if (vers == calc::v0) {
         epolar_aplus_cu<calc::V0, NON_EWALD, 1>(uind);
      } else if (vers == calc::v1) {
         epolar_aplus_cu<calc::V1, NON_EWALD, 1>(uind);
      } else if (vers == calc::v3) {
         epolar_aplus_cu<calc::V3, NON_EWALD, 1>(uind);
      } else if (vers == calc::v4) {
         epolar_aplus_cu<calc::V4, NON_EWALD, 1>(uind);
      } else if (vers == calc::v5) {
         epolar_aplus_cu<calc::V5, NON_EWALD, 1>(uind);
      } else if (vers == calc::v6) {
         epolar_aplus_cu<calc::V6, NON_EWALD, 1>(uind);
      }
   } else {
      if (vers == calc::v0) {
         epolar_aplus_cu<calc::V0, NON_EWALD, 0>(uind);
      } else if (vers == calc::v1) {
         epolar_aplus_cu<calc::V1, NON_EWALD, 0>(uind);
      } else if (vers == calc::v3) {
         epolar_aplus_cu<calc::V3, NON_EWALD, 0>(uind);
      } else if (vers == calc::v4) {
         epolar_aplus_cu<calc::V4, NON_EWALD, 0>(uind);
      } else if (vers == calc::v5) {
         epolar_aplus_cu<calc::V5, NON_EWALD, 0>(uind);
      } else if (vers == calc::v6) {
         epolar_aplus_cu<calc::V6, NON_EWALD, 0>(uind);
      }
   }
}


void epolar_aplus_ewald_real_cu(int vers, int use_cf, const real (*uind)[3])
{
   if (use_cf) {
      if (vers == calc::v0) {
         epolar_aplus_cu<calc::V0, EWALD, 1>(uind);
         //assert(false && "CFLX must compute gradient.");
      } else if (vers == calc::v1) {
         epolar_aplus_cu<calc::V1, EWALD, 1>(uind);
      } else if (vers == calc::v3) {
         epolar_aplus_cu<calc::V3, EWALD, 1>(uind);
         //assert(false && "CFLX must compute gradient.");
      } else if (vers == calc::v4) {
         epolar_aplus_cu<calc::V4, EWALD, 1>(uind);
      } else if (vers == calc::v5) {
         epolar_aplus_cu<calc::V5, EWALD, 1>(uind);
      } else if (vers == calc::v6) {
         epolar_aplus_cu<calc::V6, EWALD, 1>(uind);
      }
   } else {
      if (vers == calc::v0) {
         epolar_aplus_cu<calc::V0, EWALD, 0>(uind);
      } else if (vers == calc::v1) {
         epolar_aplus_cu<calc::V1, EWALD, 0>(uind);
      } else if (vers == calc::v3) {
         epolar_aplus_cu<calc::V3, EWALD, 0>(uind);
      } else if (vers == calc::v4) {
         epolar_aplus_cu<calc::V4, EWALD, 0>(uind);
      } else if (vers == calc::v5) {
         epolar_aplus_cu<calc::V5, EWALD, 0>(uind);
      } else if (vers == calc::v6) {
         epolar_aplus_cu<calc::V6, EWALD, 0>(uind);
      }
   }
}
}
