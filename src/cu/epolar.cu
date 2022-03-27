#include "add.h"
#include "epolar_trq.h"
#include "ff/amoeba/epolar.h"
#include "ff/energy.h"
#include "ff/image.h"
#include "ff/pme.h"
#include "ff/switch.h"
#include "launch.h"
#include "mod/elecamoeba.h"
#include "mod/elecpchg.h"
#include "mod/nblist.h"
#include "seq/pair_polar.h"
#include "seq/triangle.h"
#include "tool/gpucard.h"

namespace tinker {
// ck.py Version 2.0.2
template <class Ver, class ETYP>
__global__
void epolar_cu1(int n, TINKER_IMAGE_PARAMS, count_buffer restrict nep, energy_buffer restrict ep,
   virial_buffer restrict vep, grad_prec* restrict gx, grad_prec* restrict gy,
   grad_prec* restrict gz, real off, const unsigned* restrict mdpuinfo, int nexclude,
   const int (*restrict exclude)[2], const real (*restrict exclude_scale)[4],
   const real* restrict x, const real* restrict y, const real* restrict z,
   const Spatial::SortedAtom* restrict sorted, int nakpl, const int* restrict iakpl, int niak,
   const int* restrict iak, const int* restrict lst, real (*restrict ufld)[3],
   real (*restrict dufld)[6], const real (*restrict rpole)[10], const real (*restrict uind)[3],
   const real (*restrict uinp)[3], const real* restrict thole, const real* restrict pdamp, real f,
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
   __shared__ real uidx[BLOCK_DIM];
   __shared__ real uidy[BLOCK_DIM];
   __shared__ real uidz[BLOCK_DIM];
   __shared__ real uipx[BLOCK_DIM];
   __shared__ real uipy[BLOCK_DIM];
   __shared__ real uipz[BLOCK_DIM];
   __shared__ real pdi[BLOCK_DIM];
   __shared__ real pti[BLOCK_DIM];
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
   real ukdx;
   real ukdy;
   real ukdz;
   real ukpx;
   real ukpy;
   real ukpz;
   real pdk;
   real ptk;

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
      }

      int i = exclude[ii][0];
      int k = exclude[ii][1];
      real scaleb = exclude_scale[ii][1];
      real scalec = exclude_scale[ii][2];
      real scaled = exclude_scale[ii][3];

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
      uidx[klane] = uind[i][0];
      uidy[klane] = uind[i][1];
      uidz[klane] = uind[i][2];
      uipx[klane] = uinp[i][0];
      uipy[klane] = uinp[i][1];
      uipz[klane] = uinp[i][2];
      pdi[klane] = pdamp[i];
      pti[klane] = thole[i];
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
      ukdx = uind[k][0];
      ukdy = uind[k][1];
      ukdz = uind[k][2];
      ukpx = uinp[k][0];
      ukpy = uinp[k][1];
      ukpz = uinp[k][2];
      pdk = pdamp[k];
      ptk = thole[k];

      constexpr bool incl = true;
      real xr = xk - xi[klane];
      real yr = yk - yi[klane];
      real zr = zk - zi[klane];
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         real e, vxx, vyx, vzx, vyy, vzy, vzz;
         real e1, vxx1, vyx1, vzx1, vyy1, vzy1, vzz1;
         pair_polar_v2<Ver, ETYP>(r2, xr, yr, zr, 1, 1, 1, //
            ci[klane], dix[klane], diy[klane], diz[klane], qixx[klane], qixy[klane], qixz[klane],
            qiyy[klane], qiyz[klane], qizz[klane], uidx[klane], uidy[klane], uidz[klane],
            uipx[klane], uipy[klane], uipz[klane], pdi[klane], pti[klane], //
            ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukdx, ukdy, ukdz, ukpx, ukpy,
            ukpz, pdk, ptk, //
            f, aewald,      //
            frcxi[klane], frcyi[klane], frczi[klane], frcxk, frcyk, frczk, ufld0i[klane],
            ufld1i[klane], ufld2i[klane], ufld0k, ufld1k,
            ufld2k, //
            dufld0i[klane], dufld1i[klane], dufld2i[klane], dufld3i[klane], dufld4i[klane],
            dufld5i[klane], dufld0k, dufld1k, dufld2k, dufld3k, dufld4k, dufld5k, //
            e1, vxx1, vyx1, vzx1, vyy1, vzy1, vzz1);
         pair_polar_v2<Ver, NON_EWALD>(r2, xr, yr, zr, scaleb - 1, scalec - 1, scaled - 1, //
            ci[klane], dix[klane], diy[klane], diz[klane], qixx[klane], qixy[klane], qixz[klane],
            qiyy[klane], qiyz[klane], qizz[klane], uidx[klane], uidy[klane], uidz[klane],
            uipx[klane], uipy[klane], uipz[klane], pdi[klane], pti[klane], //
            ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukdx, ukdy, ukdz, ukpx, ukpy,
            ukpz, pdk, ptk, //
            f, aewald,      //
            frcxi[klane], frcyi[klane], frczi[klane], frcxk, frcyk, frczk, ufld0i[klane],
            ufld1i[klane], ufld2i[klane], ufld0k, ufld1k,
            ufld2k, //
            dufld0i[klane], dufld1i[klane], dufld2i[klane], dufld3i[klane], dufld4i[klane],
            dufld5i[klane], dufld0k, dufld1k, dufld2k, dufld3k, dufld4k, dufld5k, //
            e, vxx, vyx, vzx, vyy, vzy, vzz);
         if CONSTEXPR (do_e) {
            e = e + e1;
            eptl += floatTo<ebuf_prec>(e);
            if CONSTEXPR (do_a) {
               if (scalec != 0 and e != 0) // pscale != 0
                  neptl += 1;
            }
         }
         if CONSTEXPR (do_v) {
            veptlxx += floatTo<vbuf_prec>(vxx + vxx1);
            veptlyx += floatTo<vbuf_prec>(vyx + vyx1);
            veptlzx += floatTo<vbuf_prec>(vzx + vzx1);
            veptlyy += floatTo<vbuf_prec>(vyy + vyy1);
            veptlzy += floatTo<vbuf_prec>(vzy + vzy1);
            veptlzz += floatTo<vbuf_prec>(vzz + vzz1);
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
      uidx[threadIdx.x] = uind[i][0];
      uidy[threadIdx.x] = uind[i][1];
      uidz[threadIdx.x] = uind[i][2];
      uipx[threadIdx.x] = uinp[i][0];
      uipy[threadIdx.x] = uinp[i][1];
      uipz[threadIdx.x] = uinp[i][2];
      pdi[threadIdx.x] = pdamp[i];
      pti[threadIdx.x] = thole[i];
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
      ukdx = uind[k][0];
      ukdy = uind[k][1];
      ukdz = uind[k][2];
      ukpx = uinp[k][0];
      ukpy = uinp[k][1];
      ukpz = uinp[k][2];
      pdk = pdamp[k];
      ptk = thole[k];

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
            pair_polar_v2<Ver, ETYP>(r2, xr, yr, zr, 1, 1, 1, //
               ci[klane], dix[klane], diy[klane], diz[klane], qixx[klane], qixy[klane], qixz[klane],
               qiyy[klane], qiyz[klane], qizz[klane], uidx[klane], uidy[klane], uidz[klane],
               uipx[klane], uipy[klane], uipz[klane], pdi[klane], pti[klane], //
               ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukdx, ukdy, ukdz, ukpx, ukpy,
               ukpz, pdk, ptk, //
               f, aewald,      //
               frcxi[klane], frcyi[klane], frczi[klane], frcxk, frcyk, frczk, ufld0i[klane],
               ufld1i[klane], ufld2i[klane], ufld0k, ufld1k,
               ufld2k, //
               dufld0i[klane], dufld1i[klane], dufld2i[klane], dufld3i[klane], dufld4i[klane],
               dufld5i[klane], dufld0k, dufld1k, dufld2k, dufld3k, dufld4k, dufld5k, //
               e, vxx, vyx, vzx, vyy, vzy, vzz);
            if CONSTEXPR (do_e) {
               eptl += floatTo<ebuf_prec>(e);
               if CONSTEXPR (do_a) {
                  if (e != 0)
                     neptl += 1;
               }
            }
            if CONSTEXPR (do_v) {
               veptlxx += floatTo<vbuf_prec>(vxx);
               veptlyx += floatTo<vbuf_prec>(vyx);
               veptlzx += floatTo<vbuf_prec>(vzx);
               veptlyy += floatTo<vbuf_prec>(vyy);
               veptlzy += floatTo<vbuf_prec>(vzy);
               veptlzz += floatTo<vbuf_prec>(vzz);
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
         ufld1k = 0;
         ufld2k = 0;
         dufld0k = 0;
         dufld1k = 0;
         dufld2k = 0;
         dufld3k = 0;
         dufld4k = 0;
         dufld5k = 0;
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
      uidx[threadIdx.x] = uind[i][0];
      uidy[threadIdx.x] = uind[i][1];
      uidz[threadIdx.x] = uind[i][2];
      uipx[threadIdx.x] = uinp[i][0];
      uipy[threadIdx.x] = uinp[i][1];
      uipz[threadIdx.x] = uinp[i][2];
      pdi[threadIdx.x] = pdamp[i];
      pti[threadIdx.x] = thole[i];
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
      ukdx = uind[k][0];
      ukdy = uind[k][1];
      ukdz = uind[k][2];
      ukpx = uinp[k][0];
      ukpy = uinp[k][1];
      ukpz = uinp[k][2];
      pdk = pdamp[k];
      ptk = thole[k];

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
            pair_polar_v2<Ver, ETYP>(r2, xr, yr, zr, 1, 1, 1, //
               ci[klane], dix[klane], diy[klane], diz[klane], qixx[klane], qixy[klane], qixz[klane],
               qiyy[klane], qiyz[klane], qizz[klane], uidx[klane], uidy[klane], uidz[klane],
               uipx[klane], uipy[klane], uipz[klane], pdi[klane], pti[klane], //
               ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukdx, ukdy, ukdz, ukpx, ukpy,
               ukpz, pdk, ptk, //
               f, aewald,      //
               frcxi[klane], frcyi[klane], frczi[klane], frcxk, frcyk, frczk, ufld0i[klane],
               ufld1i[klane], ufld2i[klane], ufld0k, ufld1k,
               ufld2k, //
               dufld0i[klane], dufld1i[klane], dufld2i[klane], dufld3i[klane], dufld4i[klane],
               dufld5i[klane], dufld0k, dufld1k, dufld2k, dufld3k, dufld4k, dufld5k, //
               e, vxx, vyx, vzx, vyy, vzy, vzz);
            if CONSTEXPR (do_e) {
               eptl += floatTo<ebuf_prec>(e);
               if CONSTEXPR (do_a) {
                  if (e != 0)
                     neptl += 1;
               }
            }
            if CONSTEXPR (do_v) {
               veptlxx += floatTo<vbuf_prec>(vxx);
               veptlyx += floatTo<vbuf_prec>(vyx);
               veptlzx += floatTo<vbuf_prec>(vzx);
               veptlyy += floatTo<vbuf_prec>(vyy);
               veptlzy += floatTo<vbuf_prec>(vzy);
               veptlzz += floatTo<vbuf_prec>(vzz);
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
      }
   }

   if CONSTEXPR (do_a) {
      atomic_add(neptl, nep, ithread);
   }
   if CONSTEXPR (do_e) {
      atomic_add(eptl, ep, ithread);
   }
   if CONSTEXPR (do_v) {
      atomic_add(veptlxx, veptlyx, veptlzx, veptlyy, veptlzy, veptlzz, vep, ithread);
   }
}

template <class Ver, class ETYP>
void epolar_cu(const real (*uind)[3], const real (*uinp)[3])
{
   constexpr bool do_g = Ver::g;

   const auto& st = *mspatial_v2_unit;
   real off;
   if CONSTEXPR (eq<ETYP, EWALD>())
      off = switchOff(Switch::EWALD);
   else
      off = switchOff(Switch::MPOLE);

   const real f = 0.5f * electric / dielec;
   real aewald = 0;
   if CONSTEXPR (eq<ETYP, EWALD>()) {
      PMEUnit pu = ppme_unit;
      aewald = pu->aewald;
   }

   if CONSTEXPR (do_g) {
      darray::zero(g::q0, n, ufld, dufld);
   }
   int ngrid = gpuGridSize(BLOCK_DIM);
   epolar_cu1<Ver, ETYP><<<ngrid, BLOCK_DIM, 0, g::s0>>>(st.n, TINKER_IMAGE_ARGS, nep, ep, vir_ep,
      depx, depy, depz, off, st.si1.bit0, nmdpuexclude, mdpuexclude, mdpuexclude_scale, st.x, st.y,
      st.z, st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst, ufld, dufld, rpole, uind, uinp,
      thole, pdamp, f, aewald);

   // torque
   if CONSTEXPR (do_g) {
      launch_k1s(g::s0, n, epolar_trq_cu, //
         trqx, trqy, trqz, n, rpole, ufld, dufld);
   }
}

void epolar_nonewald_cu(int vers, const real (*uind)[3], const real (*uinp)[3])
{
   if (vers == calc::v0) {
      epolar_cu<calc::V0, NON_EWALD>(uind, uinp);
   } else if (vers == calc::v1) {
      epolar_cu<calc::V1, NON_EWALD>(uind, uinp);
   } else if (vers == calc::v3) {
      epolar_cu<calc::V3, NON_EWALD>(uind, uinp);
   } else if (vers == calc::v4) {
      epolar_cu<calc::V4, NON_EWALD>(uind, uinp);
   } else if (vers == calc::v5) {
      epolar_cu<calc::V5, NON_EWALD>(uind, uinp);
   } else if (vers == calc::v6) {
      epolar_cu<calc::V6, NON_EWALD>(uind, uinp);
   }
}

void epolar_ewald_real_cu(int vers, const real (*uind)[3], const real (*uinp)[3])
{
   if (vers == calc::v0) {
      epolar_cu<calc::V0, EWALD>(uind, udirp);
   } else if (vers == calc::v1) {
      epolar_cu<calc::V1, EWALD>(uind, uinp);
   } else if (vers == calc::v3) {
      epolar_cu<calc::V3, EWALD>(uind, uinp);
   } else if (vers == calc::v4) {
      epolar_cu<calc::V4, EWALD>(uind, uinp);
   } else if (vers == calc::v5) {
      epolar_cu<calc::V5, EWALD>(uind, uinp);
   } else if (vers == calc::v6) {
      epolar_cu<calc::V6, EWALD>(uind, uinp);
   }
}
}
