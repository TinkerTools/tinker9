#include "add.h"
#include "epolar.h"
#include "epolar_trq.h"
#include "glob.mplar.h"
#include "glob.spatial.h"
#include "image.h"
#include "launch.h"
#include "md.h"
#include "pme.h"
#include "seq_pair_polar.h"
#include "seq_triangle.h"
#include "switch.h"
#include "tool/gpu_card.h"


namespace tinker {
template <class Ver, class ETYP>
__global__
void epolar_cu1(int n, TINKER_IMAGE_PARAMS, count_buffer restrict nep,
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
                real (*restrict dufld)[6], const real (*restrict rpole)[10],
                const real (*restrict uind)[3], const real (*restrict uinp)[3],
                const real* restrict thole, const real* restrict pdamp, real f,
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


   __shared__ real shxi[BLOCK_DIM];
   __shared__ real shyi[BLOCK_DIM];
   __shared__ real shzi[BLOCK_DIM];
   real xk;
   real yk;
   real zk;
   __shared__ real shfrcxi[BLOCK_DIM];
   __shared__ real shfrcyi[BLOCK_DIM];
   __shared__ real shfrczi[BLOCK_DIM];
   __shared__ real shufld0i[BLOCK_DIM];
   __shared__ real shufld1i[BLOCK_DIM];
   __shared__ real shufld2i[BLOCK_DIM];
   __shared__ real shdufld0i[BLOCK_DIM];
   __shared__ real shdufld1i[BLOCK_DIM];
   __shared__ real shdufld2i[BLOCK_DIM];
   __shared__ real shdufld3i[BLOCK_DIM];
   __shared__ real shdufld4i[BLOCK_DIM];
   __shared__ real shdufld5i[BLOCK_DIM];
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
   __shared__ real shuidx[BLOCK_DIM];
   __shared__ real shuidy[BLOCK_DIM];
   __shared__ real shuidz[BLOCK_DIM];
   __shared__ real shuipx[BLOCK_DIM];
   __shared__ real shuipy[BLOCK_DIM];
   __shared__ real shuipz[BLOCK_DIM];
   __shared__ real shpdi[BLOCK_DIM];
   __shared__ real shpti[BLOCK_DIM];
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
   // exclude
   for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
      if CONSTEXPR (do_g) {
         shfrcxi[threadIdx.x] = 0;
         shfrcyi[threadIdx.x] = 0;
         shfrczi[threadIdx.x] = 0;
         shufld0i[threadIdx.x] = 0;
         shufld1i[threadIdx.x] = 0;
         shufld2i[threadIdx.x] = 0;
         shdufld0i[threadIdx.x] = 0;
         shdufld1i[threadIdx.x] = 0;
         shdufld2i[threadIdx.x] = 0;
         shdufld3i[threadIdx.x] = 0;
         shdufld4i[threadIdx.x] = 0;
         shdufld5i[threadIdx.x] = 0;
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


      int shi = exclude[ii][0];
      int k = exclude[ii][1];
      real scaleb = exclude_scale[ii][1];
      real scalec = exclude_scale[ii][2];
      real scaled = exclude_scale[ii][3];


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
      real uidx = uind[shi][0];
      real uidy = uind[shi][1];
      real uidz = uind[shi][2];
      real uipx = uinp[shi][0];
      real uipy = uinp[shi][1];
      real uipz = uinp[shi][2];
      real pdi = pdamp[shi];
      real pti = thole[shi];
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


      int klane = threadIdx.x;
      constexpr bool incl = true;
      real xr = xk - xi;
      real yr = yk - yi;
      real zr = zk - zi;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         real e, vxx, vyx, vzx, vyy, vzy, vzz;
         real e1, vxx1, vyx1, vzx1, vyy1, vzy1, vzz1;
         pair_polar_v2<Ver, ETYP>(
            r2, xr, yr, zr, 1, 1, 1, //
            ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz, uidx, uidy,
            uidz, uipx, uipy, uipz, pdi, pti, //
            ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukdx, ukdy,
            ukdz, ukpx, ukpy, ukpz, pdk, ptk, //
            f, aewald,                        //
            shfrcxi[klane], shfrcyi[klane], shfrczi[klane], frcxk, frcyk, frczk,
            shufld0i[klane], shufld1i[klane], shufld2i[klane], ufld0k, ufld1k,
            ufld2k, //
            shdufld0i[klane], shdufld1i[klane], shdufld2i[klane],
            shdufld3i[klane], shdufld4i[klane], shdufld5i[klane], dufld0k,
            dufld1k, dufld2k, dufld3k, dufld4k, dufld5k, //
            e1, vxx1, vyx1, vzx1, vyy1, vzy1, vzz1);
         pair_polar_v2<Ver, NON_EWALD>(
            r2, xr, yr, zr, scaleb - 1, scalec - 1, scaled - 1, //
            ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz, uidx, uidy,
            uidz, uipx, uipy, uipz, pdi, pti, //
            ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukdx, ukdy,
            ukdz, ukpx, ukpy, ukpz, pdk, ptk, //
            f, aewald,                        //
            shfrcxi[klane], shfrcyi[klane], shfrczi[klane], frcxk, frcyk, frczk,
            shufld0i[klane], shufld1i[klane], shufld2i[klane], ufld0k, ufld1k,
            ufld2k, //
            shdufld0i[klane], shdufld1i[klane], shdufld2i[klane],
            shdufld3i[klane], shdufld4i[klane], shdufld5i[klane], dufld0k,
            dufld1k, dufld2k, dufld3k, dufld4k, dufld5k, //
            e, vxx, vyx, vzx, vyy, vzy, vzz);
         if CONSTEXPR (do_e) {
            e = e + e1;
            eptl += cvt_to<ebuf_prec>(e);
            if CONSTEXPR (do_a) {
               if (scalec != 0 and e != 0) // pscale != 0
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
      } // end if (include)


      if CONSTEXPR (do_g) {
         atomic_add(shfrcxi[threadIdx.x], gx, shi);
         atomic_add(shfrcyi[threadIdx.x], gy, shi);
         atomic_add(shfrczi[threadIdx.x], gz, shi);
         atomic_add(shufld0i[threadIdx.x], &ufld[shi][0]);
         atomic_add(shufld1i[threadIdx.x], &ufld[shi][1]);
         atomic_add(shufld2i[threadIdx.x], &ufld[shi][2]);
         atomic_add(shdufld0i[threadIdx.x], &dufld[shi][0]);
         atomic_add(shdufld1i[threadIdx.x], &dufld[shi][1]);
         atomic_add(shdufld2i[threadIdx.x], &dufld[shi][2]);
         atomic_add(shdufld3i[threadIdx.x], &dufld[shi][3]);
         atomic_add(shdufld4i[threadIdx.x], &dufld[shi][4]);
         atomic_add(shdufld5i[threadIdx.x], &dufld[shi][5]);
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


   //* /
   // block pairs that have scale factors
   for (int iw = iwarp; iw < nakpl; iw += nwarp) {
      if CONSTEXPR (do_g) {
         shfrcxi[threadIdx.x] = 0;
         shfrcyi[threadIdx.x] = 0;
         shfrczi[threadIdx.x] = 0;
         shufld0i[threadIdx.x] = 0;
         shufld1i[threadIdx.x] = 0;
         shufld2i[threadIdx.x] = 0;
         shdufld0i[threadIdx.x] = 0;
         shdufld1i[threadIdx.x] = 0;
         shdufld2i[threadIdx.x] = 0;
         shdufld3i[threadIdx.x] = 0;
         shdufld4i[threadIdx.x] = 0;
         shdufld5i[threadIdx.x] = 0;
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
      shuidx[threadIdx.x] = uind[shi][0];
      shuidy[threadIdx.x] = uind[shi][1];
      shuidz[threadIdx.x] = uind[shi][2];
      shuipx[threadIdx.x] = uinp[shi][0];
      shuipy[threadIdx.x] = uinp[shi][1];
      shuipz[threadIdx.x] = uinp[shi][2];
      shpdi[threadIdx.x] = pdamp[shi];
      shpti[threadIdx.x] = thole[shi];
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
         real uidx = shuidx[klane];
         real uidy = shuidy[klane];
         real uidz = shuidz[klane];
         real uipx = shuipx[klane];
         real uipy = shuipy[klane];
         real uipz = shuipz[klane];
         real pdi = shpdi[klane];
         real pti = shpti[klane];


         bool incl = iid < kid and kid < n;
         incl = incl and (mdpuinfo0 & srcmask) == 0;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real e, vxx, vyx, vzx, vyy, vzy, vzz;
            pair_polar_v2<Ver, ETYP>(
               r2, xr, yr, zr, 1, 1, 1, //
               ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz, uidx,
               uidy, uidz, uipx, uipy, uipz, pdi, pti, //
               ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukdx,
               ukdy, ukdz, ukpx, ukpy, ukpz, pdk, ptk, //
               f, aewald,                              //
               shfrcxi[klane], shfrcyi[klane], shfrczi[klane], frcxk, frcyk,
               frczk, shufld0i[klane], shufld1i[klane], shufld2i[klane], ufld0k,
               ufld1k, ufld2k, //
               shdufld0i[klane], shdufld1i[klane], shdufld2i[klane],
               shdufld3i[klane], shdufld4i[klane], shdufld5i[klane], dufld0k,
               dufld1k, dufld2k, dufld3k, dufld4k, dufld5k, //
               e, vxx, vyx, vzx, vyy, vzy, vzz);
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
         } // end if (include)


         shiid = __shfl_sync(ALL_LANES, shiid, ilane + 1);
      }


      if CONSTEXPR (do_g) {
         atomic_add(shfrcxi[threadIdx.x], gx, shi);
         atomic_add(shfrcyi[threadIdx.x], gy, shi);
         atomic_add(shfrczi[threadIdx.x], gz, shi);
         atomic_add(shufld0i[threadIdx.x], &ufld[shi][0]);
         atomic_add(shufld1i[threadIdx.x], &ufld[shi][1]);
         atomic_add(shufld2i[threadIdx.x], &ufld[shi][2]);
         atomic_add(shdufld0i[threadIdx.x], &dufld[shi][0]);
         atomic_add(shdufld1i[threadIdx.x], &dufld[shi][1]);
         atomic_add(shdufld2i[threadIdx.x], &dufld[shi][2]);
         atomic_add(shdufld3i[threadIdx.x], &dufld[shi][3]);
         atomic_add(shdufld4i[threadIdx.x], &dufld[shi][4]);
         atomic_add(shdufld5i[threadIdx.x], &dufld[shi][5]);
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


   //* /
   // block-atoms
   for (int iw = iwarp; iw < niak; iw += nwarp) {
      if CONSTEXPR (do_g) {
         shfrcxi[threadIdx.x] = 0;
         shfrcyi[threadIdx.x] = 0;
         shfrczi[threadIdx.x] = 0;
         shufld0i[threadIdx.x] = 0;
         shufld1i[threadIdx.x] = 0;
         shufld2i[threadIdx.x] = 0;
         shdufld0i[threadIdx.x] = 0;
         shdufld1i[threadIdx.x] = 0;
         shdufld2i[threadIdx.x] = 0;
         shdufld3i[threadIdx.x] = 0;
         shdufld4i[threadIdx.x] = 0;
         shdufld5i[threadIdx.x] = 0;
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
      shuidx[threadIdx.x] = uind[shi][0];
      shuidy[threadIdx.x] = uind[shi][1];
      shuidz[threadIdx.x] = uind[shi][2];
      shuipx[threadIdx.x] = uinp[shi][0];
      shuipy[threadIdx.x] = uinp[shi][1];
      shuipz[threadIdx.x] = uinp[shi][2];
      shpdi[threadIdx.x] = pdamp[shi];
      shpti[threadIdx.x] = thole[shi];
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
         real uidx = shuidx[klane];
         real uidy = shuidy[klane];
         real uidz = shuidz[klane];
         real uipx = shuipx[klane];
         real uipy = shuipy[klane];
         real uipz = shuipz[klane];
         real pdi = shpdi[klane];
         real pti = shpti[klane];


         bool incl = atomk > 0;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real e, vxx, vyx, vzx, vyy, vzy, vzz;
            pair_polar_v2<Ver, ETYP>(
               r2, xr, yr, zr, 1, 1, 1, //
               ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz, uidx,
               uidy, uidz, uipx, uipy, uipz, pdi, pti, //
               ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukdx,
               ukdy, ukdz, ukpx, ukpy, ukpz, pdk, ptk, //
               f, aewald,                              //
               shfrcxi[klane], shfrcyi[klane], shfrczi[klane], frcxk, frcyk,
               frczk, shufld0i[klane], shufld1i[klane], shufld2i[klane], ufld0k,
               ufld1k, ufld2k, //
               shdufld0i[klane], shdufld1i[klane], shdufld2i[klane],
               shdufld3i[klane], shdufld4i[klane], shdufld5i[klane], dufld0k,
               dufld1k, dufld2k, dufld3k, dufld4k, dufld5k, //
               e, vxx, vyx, vzx, vyy, vzy, vzz);
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
         } // end if (include)
      }


      if CONSTEXPR (do_g) {
         atomic_add(shfrcxi[threadIdx.x], gx, shi);
         atomic_add(shfrcyi[threadIdx.x], gy, shi);
         atomic_add(shfrczi[threadIdx.x], gz, shi);
         atomic_add(shufld0i[threadIdx.x], &ufld[shi][0]);
         atomic_add(shufld1i[threadIdx.x], &ufld[shi][1]);
         atomic_add(shufld2i[threadIdx.x], &ufld[shi][2]);
         atomic_add(shdufld0i[threadIdx.x], &dufld[shi][0]);
         atomic_add(shdufld1i[threadIdx.x], &dufld[shi][1]);
         atomic_add(shdufld2i[threadIdx.x], &dufld[shi][2]);
         atomic_add(shdufld3i[threadIdx.x], &dufld[shi][3]);
         atomic_add(shdufld4i[threadIdx.x], &dufld[shi][4]);
         atomic_add(shdufld5i[threadIdx.x], &dufld[shi][5]);
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
} // generated by ComplexKernelBuilder (ck.py) 1.5.2


template <class Ver, class ETYP>
void epolar_cu(const real (*uind)[3], const real (*uinp)[3])
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
      darray::zero(PROCEED_NEW_Q, n, ufld, dufld);
   }
   int ngrid = get_grid_size(BLOCK_DIM);
   epolar_cu1<Ver, ETYP><<<ngrid, BLOCK_DIM, 0, nonblk>>>(
      st.n, TINKER_IMAGE_ARGS, nep, ep, vir_ep, depx, depy, depz, off,
      st.si1.bit0, nmdpuexclude, mdpuexclude, mdpuexclude_scale, st.x, st.y,
      st.z, st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst, ufld, dufld,
      rpole, uind, uinp, thole, pdamp, f, aewald);


   // torque
   if CONSTEXPR (do_g) {
      launch_k1s(nonblk, n, epolar_trq_cu, //
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


void epolar_ewald_real_cu(int vers, const real (*uind)[3],
                          const real (*uinp)[3])
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
