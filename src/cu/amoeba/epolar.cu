#include "ff/amoebamod.h"
#include "ff/image.h"
#include "ff/pme.h"
#include "ff/spatial.h"
#include "ff/switch.h"
#include "seq/epolartorque.h"
#include "seq/launch.h"
#include "seq/pair_polar.h"
#include "seq/triangle.h"

namespace tinker {
// ck.py Version 2.0.2
template <class Ver, class ETYP>
__global__
void epolar_cu1(int n, TINKER_IMAGE_PARAMS, CountBuffer restrict nep, EnergyBuffer restrict ep,
   VirialBuffer restrict vep, grad_prec* restrict gx, grad_prec* restrict gy,
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
   using ebuf_prec = EnergyBufferTraits::type;
   ebuf_prec eptl;
   if CONSTEXPR (do_e) {
      eptl = 0;
   }
   using vbuf_prec = VirialBufferTraits::type;
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
      ci[klane] = rpole[i][MPL_PME_0];
      dix[klane] = rpole[i][MPL_PME_X];
      diy[klane] = rpole[i][MPL_PME_Y];
      diz[klane] = rpole[i][MPL_PME_Z];
      qixx[klane] = rpole[i][MPL_PME_XX];
      qixy[klane] = rpole[i][MPL_PME_XY];
      qixz[klane] = rpole[i][MPL_PME_XZ];
      qiyy[klane] = rpole[i][MPL_PME_YY];
      qiyz[klane] = rpole[i][MPL_PME_YZ];
      qizz[klane] = rpole[i][MPL_PME_ZZ];
      uidx[klane] = uind[i][0];
      uidy[klane] = uind[i][1];
      uidz[klane] = uind[i][2];
      uipx[klane] = uinp[i][0];
      uipy[klane] = uinp[i][1];
      uipz[klane] = uinp[i][2];
      pdi[klane] = pdamp[i];
      pti[klane] = thole[i];
      ck = rpole[k][MPL_PME_0];
      dkx = rpole[k][MPL_PME_X];
      dky = rpole[k][MPL_PME_Y];
      dkz = rpole[k][MPL_PME_Z];
      qkxx = rpole[k][MPL_PME_XX];
      qkxy = rpole[k][MPL_PME_XY];
      qkxz = rpole[k][MPL_PME_XZ];
      qkyy = rpole[k][MPL_PME_YY];
      qkyz = rpole[k][MPL_PME_YZ];
      qkzz = rpole[k][MPL_PME_ZZ];
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

      ci[threadIdx.x] = rpole[i][MPL_PME_0];
      dix[threadIdx.x] = rpole[i][MPL_PME_X];
      diy[threadIdx.x] = rpole[i][MPL_PME_Y];
      diz[threadIdx.x] = rpole[i][MPL_PME_Z];
      qixx[threadIdx.x] = rpole[i][MPL_PME_XX];
      qixy[threadIdx.x] = rpole[i][MPL_PME_XY];
      qixz[threadIdx.x] = rpole[i][MPL_PME_XZ];
      qiyy[threadIdx.x] = rpole[i][MPL_PME_YY];
      qiyz[threadIdx.x] = rpole[i][MPL_PME_YZ];
      qizz[threadIdx.x] = rpole[i][MPL_PME_ZZ];
      uidx[threadIdx.x] = uind[i][0];
      uidy[threadIdx.x] = uind[i][1];
      uidz[threadIdx.x] = uind[i][2];
      uipx[threadIdx.x] = uinp[i][0];
      uipy[threadIdx.x] = uinp[i][1];
      uipz[threadIdx.x] = uinp[i][2];
      pdi[threadIdx.x] = pdamp[i];
      pti[threadIdx.x] = thole[i];
      ck = rpole[k][MPL_PME_0];
      dkx = rpole[k][MPL_PME_X];
      dky = rpole[k][MPL_PME_Y];
      dkz = rpole[k][MPL_PME_Z];
      qkxx = rpole[k][MPL_PME_XX];
      qkxy = rpole[k][MPL_PME_XY];
      qkxz = rpole[k][MPL_PME_XZ];
      qkyy = rpole[k][MPL_PME_YY];
      qkyz = rpole[k][MPL_PME_YZ];
      qkzz = rpole[k][MPL_PME_ZZ];
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

      ci[threadIdx.x] = rpole[i][MPL_PME_0];
      dix[threadIdx.x] = rpole[i][MPL_PME_X];
      diy[threadIdx.x] = rpole[i][MPL_PME_Y];
      diz[threadIdx.x] = rpole[i][MPL_PME_Z];
      qixx[threadIdx.x] = rpole[i][MPL_PME_XX];
      qixy[threadIdx.x] = rpole[i][MPL_PME_XY];
      qixz[threadIdx.x] = rpole[i][MPL_PME_XZ];
      qiyy[threadIdx.x] = rpole[i][MPL_PME_YY];
      qiyz[threadIdx.x] = rpole[i][MPL_PME_YZ];
      qizz[threadIdx.x] = rpole[i][MPL_PME_ZZ];
      uidx[threadIdx.x] = uind[i][0];
      uidy[threadIdx.x] = uind[i][1];
      uidz[threadIdx.x] = uind[i][2];
      uipx[threadIdx.x] = uinp[i][0];
      uipy[threadIdx.x] = uinp[i][1];
      uipz[threadIdx.x] = uinp[i][2];
      pdi[threadIdx.x] = pdamp[i];
      pti[threadIdx.x] = thole[i];
      ck = rpole[k][MPL_PME_0];
      dkx = rpole[k][MPL_PME_X];
      dky = rpole[k][MPL_PME_Y];
      dkz = rpole[k][MPL_PME_Z];
      qkxx = rpole[k][MPL_PME_XX];
      qkxy = rpole[k][MPL_PME_XY];
      qkxz = rpole[k][MPL_PME_XZ];
      qkyy = rpole[k][MPL_PME_YY];
      qkyz = rpole[k][MPL_PME_YZ];
      qkzz = rpole[k][MPL_PME_ZZ];
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
static void epolar_cu(const real (*uind)[3], const real (*uinp)[3])
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
      launch_k1s(g::s0, n, epolarTorque_cu, //
         trqx, trqy, trqz, n, rpole, ufld, dufld);
   }
}

void epolarNonEwald_cu(int vers, const real (*uind)[3], const real (*uinp)[3])
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

void epolarEwaldReal_cu(int vers, const real (*uind)[3], const real (*uinp)[3])
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

namespace tinker {
__global__
void epolarEwaldRecipSelfEp_cu(int n, EnergyBuffer restrict ep, real f, //
   const real (*restrict fuind)[3], const real (*fphi)[20])
{
   int ithread = ITHREAD;
   for (int i = ithread; i < n; i += STRIDE) {
      real e = 0.5f * f *
         (fuind[i][0] * fphi[i][1] + fuind[i][1] * fphi[i][2] + fuind[i][2] * fphi[i][3]);
      atomic_add(e, ep, ithread);
   }
}

__global__
void epolarEwaldRecipSelfDep_cu(int n, real f,                                   //
   grad_prec* restrict depx, grad_prec* restrict depy, grad_prec* restrict depz, //
   const real (*restrict fmp)[10], const real (*restrict fphi)[20], const real (*restrict fuind)[3],
   const real (*restrict fuinp)[3], const real (*restrict fphid)[10],
   const real (*restrict fphip)[10], const real (*restrict fphidp)[20], //
   int nfft1, int nfft2, int nfft3, TINKER_IMAGE_PARAMS)
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      // data deriv1  / 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /
      // data deriv2  / 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /
      // data deriv3  / 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /
      constexpr int deriv1[10] = {1, 4, 7, 8, 10, 15, 17, 13, 14, 19};
      constexpr int deriv2[10] = {2, 7, 5, 9, 13, 11, 18, 15, 19, 16};
      constexpr int deriv3[10] = {3, 8, 9, 6, 14, 16, 12, 19, 17, 18};

      real f1 = 0;
      real f2 = 0;
      real f3 = 0;
      for (int k = 0; k < 3; ++k) {
         int j1 = deriv1[k + 1];
         int j2 = deriv2[k + 1];
         int j3 = deriv3[k + 1];
         f1 += (fuind[i][k] + fuinp[i][k]) * fphi[i][j1];
         f2 += (fuind[i][k] + fuinp[i][k]) * fphi[i][j2];
         f3 += (fuind[i][k] + fuinp[i][k]) * fphi[i][j3];
         // if poltyp .eq. 'MUTUAL'
         f1 += fuind[i][k] * fphip[i][j1] + fuinp[i][k] * fphid[i][j1];
         f2 += fuind[i][k] * fphip[i][j2] + fuinp[i][k] * fphid[i][j2];
         f3 += fuind[i][k] * fphip[i][j3] + fuinp[i][k] * fphid[i][j3];
         // end if
      }
      for (int k = 0; k < 10; ++k) {
         f1 += fmp[i][k] * fphidp[i][deriv1[k]];
         f2 += fmp[i][k] * fphidp[i][deriv2[k]];
         f3 += fmp[i][k] * fphidp[i][deriv3[k]];
      }
      f1 *= 0.5f * nfft1;
      f2 *= 0.5f * nfft2;
      f3 *= 0.5f * nfft3;
      real h1 = recipa.x * f1 + recipb.x * f2 + recipc.x * f3;
      real h2 = recipa.y * f1 + recipb.y * f2 + recipc.y * f3;
      real h3 = recipa.z * f1 + recipb.z * f2 + recipc.z * f3;
      atomic_add(h1 * f, depx, i);
      atomic_add(h2 * f, depy, i);
      atomic_add(h3 * f, depz, i);
   }
}

template <bool do_e, bool do_g, bool do_a>
__global__
void epolarEwaldRecipSelfEptrq_cu(int n, real term, real fterm,   //
   CountBuffer restrict nep, EnergyBuffer restrict ep,            //
   real* restrict trqx, real* restrict trqy, real* restrict trqz, //
   const real (*restrict rpole)[MPL_TOTAL], const real (*restrict cmp)[10],
   const real (*restrict gpu_uind)[3], const real (*restrict gpu_uinp)[3],
   const real (*restrict cphidp)[10])
{
   int ithread = ITHREAD;
   for (int i = ithread; i < n; i += STRIDE) {
      real dix = rpole[i][MPL_PME_X];
      real diy = rpole[i][MPL_PME_Y];
      real diz = rpole[i][MPL_PME_Z];
      real uix = 0.5f * (gpu_uind[i][0] + gpu_uinp[i][0]);
      real uiy = 0.5f * (gpu_uind[i][1] + gpu_uinp[i][1]);
      real uiz = 0.5f * (gpu_uind[i][2] + gpu_uinp[i][2]);

      if CONSTEXPR (do_g) {
         real tep1 = cmp[i][3] * cphidp[i][2] - cmp[i][2] * cphidp[i][3] +
            2 * (cmp[i][6] - cmp[i][5]) * cphidp[i][9] + cmp[i][8] * cphidp[i][7] +
            cmp[i][9] * cphidp[i][5] - cmp[i][7] * cphidp[i][8] - cmp[i][9] * cphidp[i][6];
         real tep2 = cmp[i][1] * cphidp[i][3] - cmp[i][3] * cphidp[i][1] +
            2 * (cmp[i][4] - cmp[i][6]) * cphidp[i][8] + cmp[i][7] * cphidp[i][9] +
            cmp[i][8] * cphidp[i][6] - cmp[i][8] * cphidp[i][4] - cmp[i][9] * cphidp[i][7];
         real tep3 = cmp[i][2] * cphidp[i][1] - cmp[i][1] * cphidp[i][2] +
            2 * (cmp[i][5] - cmp[i][4]) * cphidp[i][7] + cmp[i][7] * cphidp[i][4] +
            cmp[i][9] * cphidp[i][8] - cmp[i][7] * cphidp[i][5] - cmp[i][8] * cphidp[i][9];

         // self term

         tep1 += term * (diy * uiz - diz * uiy);
         tep2 += term * (diz * uix - dix * uiz);
         tep3 += term * (dix * uiy - diy * uix);

         atomic_add(tep1, trqx, i);
         atomic_add(tep2, trqy, i);
         atomic_add(tep3, trqz, i);
      }

      if CONSTEXPR (do_e) {
         uix = gpu_uind[i][0];
         uiy = gpu_uind[i][1];
         uiz = gpu_uind[i][2];
         real uii = dix * uix + diy * uiy + diz * uiz;
         atomic_add(fterm * uii, ep, ithread);
      }
      if CONSTEXPR (do_a)
         atomic_add(1, nep, ithread);
   }
}

__global__
void epolarEwaldRecipSelfVirial_cu1(
   size_t size, VirialBuffer restrict vir_ep, const VirialBuffer restrict vir_m)
{
   for (size_t i = ITHREAD; i < size; i += STRIDE)
      vir_ep[0][i] -= vir_m[0][i];
}

__global__
void epolarEwaldRecipSelfVirial_cu2(int n, VirialBuffer restrict vir_ep, //
   const real (*restrict cmp)[10], const real (*restrict gpu_uind)[3],
   const real (*restrict gpu_uinp)[3], const real (*restrict fphid)[10],
   const real (*restrict fphip)[10], const real (*restrict cphi)[10],
   const real (*restrict cphidp)[10], //
   int nfft1, int nfft2, int nfft3, TINKER_IMAGE_PARAMS)
{
   int ithread = ITHREAD;
   for (int i = ithread; i < n; i += STRIDE) {
      real cphid[4], cphip[4];
      real ftc[3][3];

      // frac_to_cart

      ftc[0][0] = nfft1 * recipa.x;
      ftc[1][0] = nfft2 * recipb.x;
      ftc[2][0] = nfft3 * recipc.x;
      ftc[0][1] = nfft1 * recipa.y;
      ftc[1][1] = nfft2 * recipb.y;
      ftc[2][1] = nfft3 * recipc.y;
      ftc[0][2] = nfft1 * recipa.z;
      ftc[1][2] = nfft2 * recipb.z;
      ftc[2][2] = nfft3 * recipc.z;

      for (int j = 0; j < 3; ++j) {
         cphid[j + 1] = 0;
         cphip[j + 1] = 0;
         for (int k = 0; k < 3; ++k) {
            cphid[j + 1] += ftc[k][j] * fphid[i][k + 1];
            cphip[j + 1] += ftc[k][j] * fphip[i][k + 1];
         }
      }

      real vxx = 0;
      real vyy = 0;
      real vzz = 0;
      real vxy = 0;
      real vxz = 0;
      real vyz = 0;

      vxx =
         vxx - cmp[i][1] * cphidp[i][1] - 0.5f * ((gpu_uind[i][0] + gpu_uinp[i][0]) * cphi[i][1]);
      vxy = vxy - 0.5f * (cphidp[i][1] * cmp[i][2] + cphidp[i][2] * cmp[i][1]) -
         0.25f *
            ((gpu_uind[i][1] + gpu_uinp[i][1]) * cphi[i][1] +
               (gpu_uind[i][0] + gpu_uinp[i][0]) * cphi[i][2]);
      vxz = vxz - 0.5f * (cphidp[i][1] * cmp[i][3] + cphidp[i][3] * cmp[i][1]) -
         0.25f *
            ((gpu_uind[i][2] + gpu_uinp[i][2]) * cphi[i][1] +
               (gpu_uind[i][0] + gpu_uinp[i][0]) * cphi[i][3]);
      vyy =
         vyy - cmp[i][2] * cphidp[i][2] - 0.5f * ((gpu_uind[i][1] + gpu_uinp[i][1]) * cphi[i][2]);
      vyz = vyz - 0.5f * (cphidp[i][2] * cmp[i][3] + cphidp[i][3] * cmp[i][2]) -
         0.25f *
            ((gpu_uind[i][2] + gpu_uinp[i][2]) * cphi[i][2] +
               (gpu_uind[i][1] + gpu_uinp[i][1]) * cphi[i][3]);
      vzz =
         vzz - cmp[i][3] * cphidp[i][3] - 0.5f * ((gpu_uind[i][2] + gpu_uinp[i][2]) * cphi[i][3]);
      vxx =
         vxx - 2 * cmp[i][4] * cphidp[i][4] - cmp[i][7] * cphidp[i][7] - cmp[i][8] * cphidp[i][8];
      vxy = vxy - (cmp[i][4] + cmp[i][5]) * cphidp[i][7] -
         0.5f *
            (cmp[i][7] * (cphidp[i][5] + cphidp[i][4]) + cmp[i][8] * cphidp[i][9] +
               cmp[i][9] * cphidp[i][8]);
      vxz = vxz - (cmp[i][4] + cmp[i][6]) * cphidp[i][8] -
         0.5f *
            (cmp[i][8] * (cphidp[i][4] + cphidp[i][6]) + cmp[i][7] * cphidp[i][9] +
               cmp[i][9] * cphidp[i][7]);
      vyy =
         vyy - 2 * cmp[i][5] * cphidp[i][5] - cmp[i][7] * cphidp[i][7] - cmp[i][9] * cphidp[i][9];
      vyz = vyz - (cmp[i][5] + cmp[i][6]) * cphidp[i][9] -
         0.5f *
            (cmp[i][9] * (cphidp[i][5] + cphidp[i][6]) + cmp[i][7] * cphidp[i][8] +
               cmp[i][8] * cphidp[i][7]);
      vzz =
         vzz - 2 * cmp[i][6] * cphidp[i][6] - cmp[i][8] * cphidp[i][8] - cmp[i][9] * cphidp[i][9];

      // if (poltyp .eq. 'MUTUAL')
      vxx = vxx - 0.5f * (cphid[1] * gpu_uinp[i][0] + cphip[1] * gpu_uind[i][0]);
      vxy = vxy -
         0.25f *
            (cphid[1] * gpu_uinp[i][1] + cphip[1] * gpu_uind[i][1] + cphid[2] * gpu_uinp[i][0] +
               cphip[2] * gpu_uind[i][0]);
      vxz = vxz -
         0.25f *
            (cphid[1] * gpu_uinp[i][2] + cphip[1] * gpu_uind[i][2] + cphid[3] * gpu_uinp[i][0] +
               cphip[3] * gpu_uind[i][0]);
      vyy = vyy - 0.5f * (cphid[2] * gpu_uinp[i][1] + cphip[2] * gpu_uind[i][1]);
      vyz = vyz -
         0.25f *
            (cphid[2] * gpu_uinp[i][2] + cphip[2] * gpu_uind[i][2] + cphid[3] * gpu_uinp[i][1] +
               cphip[3] * gpu_uind[i][1]);
      vzz = vzz - 0.5f * (cphid[3] * gpu_uinp[i][2] + cphip[3] * gpu_uind[i][2]);
      // end if

      atomic_add(vxx, vxy, vxz, vyy, vyz, vzz, vir_ep, ithread);
   }
}

__global__
void epolarEwaldRecipSelfVirial_cu3(
   int n, real (*restrict cmp)[10], const real (*restrict gpu_uinp)[3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      cmp[i][1] += gpu_uinp[i][0];
      cmp[i][2] += gpu_uinp[i][1];
      cmp[i][3] += gpu_uinp[i][2];
   }
}

__global__
void epolarEwaldRecipSelfVirial_cu4(int n, real (*restrict cmp)[10],
   const real (*restrict gpu_uind)[3], const real (*restrict gpu_uinp)[3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      cmp[i][1] += (gpu_uind[i][0] - gpu_uinp[i][0]);
      cmp[i][2] += (gpu_uind[i][1] - gpu_uinp[i][1]);
      cmp[i][3] += (gpu_uind[i][2] - gpu_uinp[i][2]);
   }
}

__global__
void epolarEwaldRecipSelfVirial_cu5(int ntot, int nff, VirialBuffer restrict vir_ep, //
   real f, real volterm, real pterm, const PME* restrict d, const PME* restrict p,   //
   int nfft1, int nfft2, int nfft3, TINKER_IMAGE_PARAMS)
{
   int ithread = ITHREAD;
   for (int i = ithread; i < ntot; i += STRIDE) {
      if (i == 0) {
         continue;
      }

      // const real volterm = pi * box_volume;

      int k3 = i / nff;
      int j = i - k3 * nff;
      int k2 = j / nfft1;
      int k1 = j - k2 * nfft1;

      int r1 = (k1 < (nfft1 + 1) / 2) ? k1 : (k1 - nfft1);
      int r2 = (k2 < (nfft2 + 1) / 2) ? k2 : (k2 - nfft2);
      int r3 = (k3 < (nfft3 + 1) / 2) ? k3 : (k3 - nfft3);

      real h1 = recipa.x * r1 + recipb.x * r2 + recipc.x * r3;
      real h2 = recipa.y * r1 + recipb.y * r2 + recipc.y * r3;
      real h3 = recipa.z * r1 + recipb.z * r2 + recipc.z * r3;
      real hsq = h1 * h1 + h2 * h2 + h3 * h3;
      real term = -pterm * hsq;
      real expterm = 0;
      if (term > -50) {
         real denom = volterm * hsq * d->bsmod1[k1] * d->bsmod2[k2] * d->bsmod3[k3];
         expterm = REAL_EXP(term) / denom;
         if (box_shape == BoxShape::UNBOUND)
            expterm *= (1 - REAL_COS(pi * lvec1.x * REAL_SQRT(hsq)));
         else if (box_shape == BoxShape::OCT)
            if ((k1 + k2 + k3) & 1)
               expterm = 0; // end if ((k1 + k2 + k3) % 2 != 0)

         real struc2 =
            d->qgrid[2 * i] * p->qgrid[2 * i] + d->qgrid[2 * i + 1] * p->qgrid[2 * i + 1];
         real eterm = 0.5f * f * expterm * struc2;
         real vterm = (2 / hsq) * (1 - term) * eterm;

         real vxx = (h1 * h1 * vterm - eterm);
         real vxy = h1 * h2 * vterm;
         real vxz = h1 * h3 * vterm;
         real vyy = (h2 * h2 * vterm - eterm);
         real vyz = h2 * h3 * vterm;
         real vzz = (h3 * h3 * vterm - eterm);

         atomic_add(vxx, vxy, vxz, vyy, vyz, vzz, vir_ep, ithread);
      }
   }
}

template <class Ver>
static void epolarEwaldRecipSelf_cu1(const real (*gpu_uind)[3], const real (*gpu_uinp)[3])
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;

   const PMEUnit pu = ppme_unit;
   const auto& st = *pu;
   const int nfft1 = st.nfft1;
   const int nfft2 = st.nfft2;
   const int nfft3 = st.nfft3;
   const real aewald = st.aewald;

   const real f = electric / dielec;

   auto* fphid = fdip_phi1;
   auto* fphip = fdip_phi2;

   cuindToFuind(pu, gpu_uind, gpu_uinp, fuind, fuinp);
   if CONSTEXPR (do_e) {
      launch_k1b(g::s0, n, epolarEwaldRecipSelfEp_cu, n, ep, f, fuind, fphi);
   }
   gridUind(pu, fuind, fuinp);
   fftfront(pu);
   // TODO: store vs. recompute qfac
   pmeConv(pu);
   fftback(pu);
   fphiUind(pu, fphid, fphip, fphidp);

   // increment the dipole polarization gradient contributions

   if CONSTEXPR (do_g) {
      launch_k1s(g::s0, n, epolarEwaldRecipSelfDep_cu,  //
         n, f, depx, depy, depz,                        //
         fmp, fphi, fuind, fuinp, fphid, fphip, fphidp, //
         nfft1, nfft2, nfft3, TINKER_IMAGE_ARGS);
   }

   // set the potential to be the induced dipole average

   // see also subroutine eprecip1 in epolar1.f
   // do i = 1, npole
   //    do j = 1, 10
   //       fphidp(j,i) = 0.5d0 * fphidp(j,i)
   //    end do
   // end do
   // Notice that only 10 * n elements were scaled in the original code.
   darray::scale(g::q0, n, 0.5f * f, fphidp);
   fphiToCphi(pu, fphidp, cphidp);

   // recip and self torques

   real term = f * aewald * aewald * aewald * 4 / 3 / sqrtpi;
   real fterm = -2 * f * aewald * aewald * aewald / 3 / sqrtpi;
   launch_k1b(g::s0, n, epolarEwaldRecipSelfEptrq_cu<do_e, do_g, do_a>, //
      n, term, fterm, nep, ep, trqx, trqy, trqz,                        //
      rpole, cmp, gpu_uind, gpu_uinp, cphidp);

   // recip virial

   if CONSTEXPR (do_v) {
      auto size = bufferSize() * VirialBufferTraits::value;
      launch_k1s(g::s0, size, epolarEwaldRecipSelfVirial_cu1, size, vir_ep, vir_m);

      darray::scale(g::q0, n, f, cphi, fphid, fphip);

      launch_k1b(g::s0, n, epolarEwaldRecipSelfVirial_cu2,    //
         n, vir_ep,                                           //
         cmp, gpu_uind, gpu_uinp, fphid, fphip, cphi, cphidp, //
         nfft1, nfft2, nfft3, TINKER_IMAGE_ARGS);

      // qgrip: pvu_qgrid
      const PMEUnit pvu = pvpme_unit;
      launch_k1s(g::s0, n, epolarEwaldRecipSelfVirial_cu3, n, cmp, gpu_uinp);
      cmpToFmp(pvu, cmp, fmp);
      gridMpole(pvu, fmp);
      fftfront(pvu);

      // qgrid: pu_qgrid
      launch_k1s(g::s0, n, epolarEwaldRecipSelfVirial_cu4, n, cmp, gpu_uind, gpu_uinp);
      cmpToFmp(pu, cmp, fmp);
      gridMpole(pu, fmp);
      fftfront(pu);

      const auto* d = pu.deviceptr();
      const auto* p = pvu.deviceptr();
      const int nff = nfft1 * nfft2;
      const int ntot = nfft1 * nfft2 * nfft3;
      real pterm = (pi / aewald) * (pi / aewald);
      real box_volume = boxVolume();
      real volterm = pi * box_volume;
      launch_k1b(g::s0, ntot, epolarEwaldRecipSelfVirial_cu5, //
         ntot, nff, vir_ep, f, volterm, pterm, d, p,          //
         nfft1, nfft2, nfft3, TINKER_IMAGE_ARGS);
   }
}

void epolarEwaldRecipSelf_cu(int vers, const real (*uind)[3], const real (*uinp)[3])
{
   if (vers == calc::v0) {
      epolarEwaldRecipSelf_cu1<calc::V0>(uind, uinp);
   } else if (vers == calc::v1) {
      epolarEwaldRecipSelf_cu1<calc::V1>(uind, uinp);
   } else if (vers == calc::v3) {
      epolarEwaldRecipSelf_cu1<calc::V3>(uind, uinp);
   } else if (vers == calc::v4) {
      epolarEwaldRecipSelf_cu1<calc::V4>(uind, uinp);
   } else if (vers == calc::v5) {
      epolarEwaldRecipSelf_cu1<calc::V5>(uind, uinp);
   } else if (vers == calc::v6) {
      epolarEwaldRecipSelf_cu1<calc::V6>(uind, uinp);
   }
}
}
