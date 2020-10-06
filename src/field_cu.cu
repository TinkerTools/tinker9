#include "add.h"
#include "epolar.h"
#include "glob.mplar.h"
#include "glob.spatial.h"
#include "image.h"
#include "launch.h"
#include "md.h"
#include "pme.h"
#include "seq_pair_field.h"
#include "seq_triangle.h"
#include "switch.h"
#include "tool/cudalib.h"
#include "tool/gpu_card.h"


namespace tinker {
#define DFIELDPARAS                                                            \
   real(*restrict field)[3], real(*restrict fieldp)[3],                        \
      const real *restrict thole, const real *restrict pdamp,                  \
      const real(*restrict rpole)[10], TINKER_IMAGE_PARAMS, real off2


template <class ETYP>
__launch_bounds__(BLOCK_DIM) __global__
void dfield_cu1(DFIELDPARAS, int n, const Spatial::SortedAtom* restrict sorted,
                int niak, const int* restrict iak, const int* restrict lst,
                real aewald)
{
   const int iwarp = (threadIdx.x + blockIdx.x * blockDim.x) / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);


   struct Data
   {
      real3 fkd, fkp;
      real3 rk;
      real ck;
      real dkx, dky, dkz;
      real qkxx, qkxy, qkxz, qkyy, qkyz, qkzz;
      real pdk, ptk;
   };
   __shared__ Data data[BLOCK_DIM];


   for (int iw = iwarp; iw < niak; iw += nwarp) {
      real3 fid = make_real3(0, 0, 0);
      real3 fip = make_real3(0, 0, 0);
      int atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
      real3 ri = make_real3(sorted[atomi].x, sorted[atomi].y, sorted[atomi].z);
      int i = sorted[atomi].unsorted;
      real ci = rpole[i][mpl_pme_0];
      real dix = rpole[i][mpl_pme_x];
      real diy = rpole[i][mpl_pme_y];
      real diz = rpole[i][mpl_pme_z];
      real qixx = rpole[i][mpl_pme_xx];
      real qixy = rpole[i][mpl_pme_xy];
      real qixz = rpole[i][mpl_pme_xz];
      real qiyy = rpole[i][mpl_pme_yy];
      real qiyz = rpole[i][mpl_pme_yz];
      real qizz = rpole[i][mpl_pme_zz];
      real pdi = pdamp[i];
      real pti = thole[i];


      data[threadIdx.x].fkd = make_real3(0, 0, 0);
      data[threadIdx.x].fkp = make_real3(0, 0, 0);
      int shatomk = lst[iw * WARP_SIZE + ilane];
      data[threadIdx.x].rk =
         make_real3(sorted[shatomk].x, sorted[shatomk].y, sorted[shatomk].z);
      int shk = sorted[shatomk].unsorted;
      data[threadIdx.x].ck = rpole[shk][mpl_pme_0];
      data[threadIdx.x].dkx = rpole[shk][mpl_pme_x];
      data[threadIdx.x].dky = rpole[shk][mpl_pme_y];
      data[threadIdx.x].dkz = rpole[shk][mpl_pme_z];
      data[threadIdx.x].qkxx = rpole[shk][mpl_pme_xx];
      data[threadIdx.x].qkxy = rpole[shk][mpl_pme_xy];
      data[threadIdx.x].qkxz = rpole[shk][mpl_pme_xz];
      data[threadIdx.x].qkyy = rpole[shk][mpl_pme_yy];
      data[threadIdx.x].qkyz = rpole[shk][mpl_pme_yz];
      data[threadIdx.x].qkzz = rpole[shk][mpl_pme_zz];
      data[threadIdx.x].pdk = pdamp[shk];
      data[threadIdx.x].ptk = thole[shk];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
         real3 dr = data[klane].rk - ri;


         real r2 = image2(dr.x, dr.y, dr.z);
         if (atomi < atomk && r2 <= off2) {
            if CONSTEXPR (eq<ETYP, EWALD>()) {
               pair_dfield<EWALD>(
                  r2, dr.x, dr.y, dr.z, 1, 1, ci, dix, diy, diz, qixx, qixy,
                  qixz, qiyy, qiyz, qizz, pdi, pti, data[klane].ck,
                  data[klane].dkx, data[klane].dky, data[klane].dkz,
                  data[klane].qkxx, data[klane].qkxy, data[klane].qkxz,
                  data[klane].qkyy, data[klane].qkyz, data[klane].qkzz,
                  data[klane].pdk, data[klane].ptk, aewald, fid, fip,
                  data[klane].fkd, data[klane].fkp);
            }
            if CONSTEXPR (eq<ETYP, NON_EWALD>()) {
               pair_dfield<NON_EWALD>(
                  r2, dr.x, dr.y, dr.z, 1, 1, ci, dix, diy, diz, qixx, qixy,
                  qixz, qiyy, qiyz, qizz, pdi, pti, data[klane].ck,
                  data[klane].dkx, data[klane].dky, data[klane].dkz,
                  data[klane].qkxx, data[klane].qkxy, data[klane].qkxz,
                  data[klane].qkyy, data[klane].qkyz, data[klane].qkzz,
                  data[klane].pdk, data[klane].ptk, 0, fid, fip,
                  data[klane].fkd, data[klane].fkp);
            }
         } // end if (include)
      }


      atomic_add(fid.x, &field[i][0]);
      atomic_add(fid.y, &field[i][1]);
      atomic_add(fid.z, &field[i][2]);
      atomic_add(fip.x, &fieldp[i][0]);
      atomic_add(fip.y, &fieldp[i][1]);
      atomic_add(fip.z, &fieldp[i][2]);
      atomic_add(data[threadIdx.x].fkd.x, &field[shk][0]);
      atomic_add(data[threadIdx.x].fkd.y, &field[shk][1]);
      atomic_add(data[threadIdx.x].fkd.z, &field[shk][2]);
      atomic_add(data[threadIdx.x].fkp.x, &fieldp[shk][0]);
      atomic_add(data[threadIdx.x].fkp.y, &fieldp[shk][1]);
      atomic_add(data[threadIdx.x].fkp.z, &fieldp[shk][2]);
   } // end for (iw)
}


__global__
void dfield_cu2(DFIELDPARAS, const real* restrict x, const real* restrict y,
                const real* restrict z, int ndpexclude,
                const int (*restrict dpexclude)[2],
                const real (*restrict dpexclude_scale)[2])
{
   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < ndpexclude;
        ii += blockDim.x * gridDim.x) {
      int i = dpexclude[ii][0];
      int k = dpexclude[ii][1];
      real dscale = dpexclude_scale[ii][0] - 1;
      real pscale = dpexclude_scale[ii][1] - 1;


      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real ci = rpole[i][mpl_pme_0];
      real dix = rpole[i][mpl_pme_x];
      real diy = rpole[i][mpl_pme_y];
      real diz = rpole[i][mpl_pme_z];
      real qixx = rpole[i][mpl_pme_xx];
      real qixy = rpole[i][mpl_pme_xy];
      real qixz = rpole[i][mpl_pme_xz];
      real qiyy = rpole[i][mpl_pme_yy];
      real qiyz = rpole[i][mpl_pme_yz];
      real qizz = rpole[i][mpl_pme_zz];
      real pdi = pdamp[i];
      real pti = thole[i];


      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;


      real r2 = image2(xr, yr, zr);
      if (r2 <= off2) {
         real3 fid = make_real3(0, 0, 0);
         real3 fip = make_real3(0, 0, 0);
         real3 fkd = make_real3(0, 0, 0);
         real3 fkp = make_real3(0, 0, 0);
         pair_dfield<NON_EWALD>(
            r2, xr, yr, zr, dscale, pscale, ci, dix, diy, diz, qixx, qixy, qixz,
            qiyy, qiyz, qizz, pdi, pti, rpole[k][mpl_pme_0],
            rpole[k][mpl_pme_x], rpole[k][mpl_pme_y], rpole[k][mpl_pme_z],
            rpole[k][mpl_pme_xx], rpole[k][mpl_pme_xy], rpole[k][mpl_pme_xz],
            rpole[k][mpl_pme_yy], rpole[k][mpl_pme_yz], rpole[k][mpl_pme_zz],
            pdamp[k], thole[k], 0, fid, fip, fkd, fkp);


         atomic_add(fid.x, &field[i][0]);
         atomic_add(fid.y, &field[i][1]);
         atomic_add(fid.z, &field[i][2]);
         atomic_add(fip.x, &fieldp[i][0]);
         atomic_add(fip.y, &fieldp[i][1]);
         atomic_add(fip.z, &fieldp[i][2]);


         atomic_add(fkd.x, &field[k][0]);
         atomic_add(fkd.y, &field[k][1]);
         atomic_add(fkd.z, &field[k][2]);
         atomic_add(fkp.x, &fieldp[k][0]);
         atomic_add(fkp.y, &fieldp[k][1]);
         atomic_add(fkp.z, &fieldp[k][2]);
      } // end if (include)
   }
}


void dfield_ewald_real_cu(real (*field)[3], real (*fieldp)[3])
{
   const auto& st = *mspatial_unit;
   const real off = switch_off(switch_ewald);
   const real off2 = off * off;


   const PMEUnit pu = ppme_unit;
   const real aewald = pu->aewald;
   if (st.niak > 0) {
      launch_k1s(nonblk, WARP_SIZE * st.niak, dfield_cu1<EWALD>, field, fieldp,
                 thole, pdamp, rpole, TINKER_IMAGE_ARGS, off2, n, st.sorted,
                 st.niak, st.iak, st.lst, aewald);
   }
   if (ndpexclude > 0) {
      launch_k1s(nonblk, ndpexclude, dfield_cu2, field, fieldp, thole, pdamp,
                 rpole, TINKER_IMAGE_ARGS, off2, x, y, z, ndpexclude, dpexclude,
                 dpexclude_scale);
   }
}


void dfield_nonewald_cu(real (*field)[3], real (*fieldp)[3])
{
   const auto& st = *mspatial_unit;
   const real off = switch_off(switch_mpole);
   const real off2 = off * off;


   darray::zero(PROCEED_NEW_Q, n, field, fieldp);
   if (st.niak > 0) {
      launch_k1s(nonblk, WARP_SIZE * st.niak, dfield_cu1<NON_EWALD>, field,
                 fieldp, thole, pdamp, rpole, TINKER_IMAGE_ARGS, off2, n,
                 st.sorted, st.niak, st.iak, st.lst, 0);
   }
   if (ndpexclude > 0) {
      launch_k1s(nonblk, ndpexclude, dfield_cu2, field, fieldp, thole, pdamp,
                 rpole, TINKER_IMAGE_ARGS, off2, x, y, z, ndpexclude, dpexclude,
                 dpexclude_scale);
   }
}


template <class ETYP>
__launch_bounds__(BLOCK_DIM) __global__
void ufield_cu1(int n, TINKER_IMAGE_PARAMS, real off,
                const unsigned* restrict minfo, const unsigned* restrict dpinfo,
                const unsigned* restrict uinfo, int nexclude,
                const int (*restrict exclude)[2],
                const real (*restrict exclude_scale)[4], const real* restrict x,
                const real* restrict y, const real* restrict z,
                const Spatial::SortedAtom* restrict sorted, int nakpl,
                const int* restrict iakpl, int niak, const int* restrict iak,
                const int* restrict lst, const real (*restrict uind)[3],
                const real (*restrict uinp)[3], real (*restrict field)[3],
                real (*restrict fieldp)[3], const real* restrict thole,
                const real* restrict pdamp, real aewald)
{
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);


   __shared__ real shxi[BLOCK_DIM];
   __shared__ real shyi[BLOCK_DIM];
   __shared__ real shzi[BLOCK_DIM];
   real xk;
   real yk;
   real zk;
   __shared__ real shfidx[BLOCK_DIM];
   __shared__ real shfidy[BLOCK_DIM];
   __shared__ real shfidz[BLOCK_DIM];
   __shared__ real shfipx[BLOCK_DIM];
   __shared__ real shfipy[BLOCK_DIM];
   __shared__ real shfipz[BLOCK_DIM];
   real fkdx;
   real fkdy;
   real fkdz;
   real fkpx;
   real fkpy;
   real fkpz;
   __shared__ real shuidx[BLOCK_DIM];
   __shared__ real shuidy[BLOCK_DIM];
   __shared__ real shuidz[BLOCK_DIM];
   __shared__ real shuipx[BLOCK_DIM];
   __shared__ real shuipy[BLOCK_DIM];
   __shared__ real shuipz[BLOCK_DIM];
   __shared__ real shpdi[BLOCK_DIM];
   __shared__ real shpti[BLOCK_DIM];
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
      shfidx[threadIdx.x] = 0;
      shfidy[threadIdx.x] = 0;
      shfidz[threadIdx.x] = 0;
      shfipx[threadIdx.x] = 0;
      shfipy[threadIdx.x] = 0;
      shfipz[threadIdx.x] = 0;
      fkdx = 0;
      fkdy = 0;
      fkdz = 0;
      fkpx = 0;
      fkpy = 0;
      fkpz = 0;


      int shi = exclude[ii][0];
      int k = exclude[ii][1];
      real scaled = exclude_scale[ii][3];


      real xi = x[shi];
      real yi = y[shi];
      real zi = z[shi];
      xk = x[k];
      yk = y[k];
      zk = z[k];
      real uidx = uind[shi][0];
      real uidy = uind[shi][1];
      real uidz = uind[shi][2];
      real uipx = uinp[shi][0];
      real uipy = uinp[shi][1];
      real uipz = uinp[shi][2];
      real pdi = pdamp[shi];
      real pti = thole[shi];
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
         pair_ufield_v2<ETYP>(
            r2, xr, yr, zr, scaled, aewald, uidx, uidy, uidz, uipx, uipy, uipz,
            pdi, pti, ukdx, ukdy, ukdz, ukpx, ukpy, ukpz, pdk, ptk,
            shfidx[klane], shfidy[klane], shfidz[klane], shfipx[klane],
            shfipy[klane], shfipz[klane], fkdx, fkdy, fkdz, fkpx, fkpy, fkpz);
      } // end if (include)


      atomic_add(shfidx[threadIdx.x], &field[shi][0]);
      atomic_add(shfidy[threadIdx.x], &field[shi][1]);
      atomic_add(shfidz[threadIdx.x], &field[shi][2]);
      atomic_add(shfipx[threadIdx.x], &fieldp[shi][0]);
      atomic_add(shfipy[threadIdx.x], &fieldp[shi][1]);
      atomic_add(shfipz[threadIdx.x], &fieldp[shi][2]);
      atomic_add(fkdx, &field[k][0]);
      atomic_add(fkdy, &field[k][1]);
      atomic_add(fkdz, &field[k][2]);
      atomic_add(fkpx, &fieldp[k][0]);
      atomic_add(fkpy, &fieldp[k][1]);
      atomic_add(fkpz, &fieldp[k][2]);
   }
   // */


   //* /
   // block pairs that have scale factors
   for (int iw = iwarp; iw < nakpl; iw += nwarp) {
      shfidx[threadIdx.x] = 0;
      shfidy[threadIdx.x] = 0;
      shfidz[threadIdx.x] = 0;
      shfipx[threadIdx.x] = 0;
      shfipy[threadIdx.x] = 0;
      shfipz[threadIdx.x] = 0;
      fkdx = 0;
      fkdy = 0;
      fkdz = 0;
      fkpx = 0;
      fkpy = 0;
      fkpz = 0;


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


      shuidx[threadIdx.x] = uind[shi][0];
      shuidy[threadIdx.x] = uind[shi][1];
      shuidz[threadIdx.x] = uind[shi][2];
      shuipx[threadIdx.x] = uinp[shi][0];
      shuipy[threadIdx.x] = uinp[shi][1];
      shuipz[threadIdx.x] = uinp[shi][2];
      shpdi[threadIdx.x] = pdamp[shi];
      shpti[threadIdx.x] = thole[shi];
      ukdx = uind[k][0];
      ukdy = uind[k][1];
      ukdz = uind[k][2];
      ukpx = uinp[k][0];
      ukpy = uinp[k][1];
      ukpz = uinp[k][2];
      pdk = pdamp[k];
      ptk = thole[k];


      unsigned int minfo0 = minfo[iw * WARP_SIZE + ilane];
      unsigned int dpinfo0 = dpinfo[iw * WARP_SIZE + ilane];
      unsigned int uinfo0 = uinfo[iw * WARP_SIZE + ilane];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int srcmask = 1 << srclane;
         int klane = srclane + threadIdx.x - ilane;
         int iid = shiid;
         real xi = shxi[klane];
         real yi = shyi[klane];
         real zi = shzi[klane];
         real uidx = shuidx[klane];
         real uidy = shuidy[klane];
         real uidz = shuidz[klane];
         real uipx = shuipx[klane];
         real uipy = shuipy[klane];
         real uipz = shuipz[klane];
         real pdi = shpdi[klane];
         real pti = shpti[klane];


         bool incl = iid < kid and kid < n;
         incl = incl and (minfo0 & srcmask) == 0 and
            (dpinfo0 & srcmask) == 0 and (uinfo0 & srcmask) == 0;
         real scaled = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_ufield_v2<ETYP>(r2, xr, yr, zr, scaled, aewald, uidx, uidy,
                                 uidz, uipx, uipy, uipz, pdi, pti, ukdx, ukdy,
                                 ukdz, ukpx, ukpy, ukpz, pdk, ptk,
                                 shfidx[klane], shfidy[klane], shfidz[klane],
                                 shfipx[klane], shfipy[klane], shfipz[klane],
                                 fkdx, fkdy, fkdz, fkpx, fkpy, fkpz);
         } // end if (include)


         shiid = __shfl_sync(ALL_LANES, shiid, ilane + 1);
      }


      atomic_add(shfidx[threadIdx.x], &field[shi][0]);
      atomic_add(shfidy[threadIdx.x], &field[shi][1]);
      atomic_add(shfidz[threadIdx.x], &field[shi][2]);
      atomic_add(shfipx[threadIdx.x], &fieldp[shi][0]);
      atomic_add(shfipy[threadIdx.x], &fieldp[shi][1]);
      atomic_add(shfipz[threadIdx.x], &fieldp[shi][2]);
      atomic_add(fkdx, &field[k][0]);
      atomic_add(fkdy, &field[k][1]);
      atomic_add(fkdz, &field[k][2]);
      atomic_add(fkpx, &fieldp[k][0]);
      atomic_add(fkpy, &fieldp[k][1]);
      atomic_add(fkpz, &fieldp[k][2]);
   }
   // */


   //* /
   // block-atoms
   for (int iw = iwarp; iw < niak; iw += nwarp) {
      shfidx[threadIdx.x] = 0;
      shfidy[threadIdx.x] = 0;
      shfidz[threadIdx.x] = 0;
      shfipx[threadIdx.x] = 0;
      shfipy[threadIdx.x] = 0;
      shfipz[threadIdx.x] = 0;
      fkdx = 0;
      fkdy = 0;
      fkdz = 0;
      fkpx = 0;
      fkpy = 0;
      fkpz = 0;


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


      shuidx[threadIdx.x] = uind[shi][0];
      shuidy[threadIdx.x] = uind[shi][1];
      shuidz[threadIdx.x] = uind[shi][2];
      shuipx[threadIdx.x] = uinp[shi][0];
      shuipy[threadIdx.x] = uinp[shi][1];
      shuipz[threadIdx.x] = uinp[shi][2];
      shpdi[threadIdx.x] = pdamp[shi];
      shpti[threadIdx.x] = thole[shi];
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
         real uidx = shuidx[klane];
         real uidy = shuidy[klane];
         real uidz = shuidz[klane];
         real uipx = shuipx[klane];
         real uipy = shuipy[klane];
         real uipz = shuipz[klane];
         real pdi = shpdi[klane];
         real pti = shpti[klane];


         bool incl = atomk > 0;
         real scaled = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_ufield_v2<ETYP>(r2, xr, yr, zr, scaled, aewald, uidx, uidy,
                                 uidz, uipx, uipy, uipz, pdi, pti, ukdx, ukdy,
                                 ukdz, ukpx, ukpy, ukpz, pdk, ptk,
                                 shfidx[klane], shfidy[klane], shfidz[klane],
                                 shfipx[klane], shfipy[klane], shfipz[klane],
                                 fkdx, fkdy, fkdz, fkpx, fkpy, fkpz);
         } // end if (include)
      }


      atomic_add(shfidx[threadIdx.x], &field[shi][0]);
      atomic_add(shfidy[threadIdx.x], &field[shi][1]);
      atomic_add(shfidz[threadIdx.x], &field[shi][2]);
      atomic_add(shfipx[threadIdx.x], &fieldp[shi][0]);
      atomic_add(shfipy[threadIdx.x], &fieldp[shi][1]);
      atomic_add(shfipz[threadIdx.x], &fieldp[shi][2]);
      atomic_add(fkdx, &field[k][0]);
      atomic_add(fkdy, &field[k][1]);
      atomic_add(fkdz, &field[k][2]);
      atomic_add(fkpx, &fieldp[k][0]);
      atomic_add(fkpy, &fieldp[k][1]);
      atomic_add(fkpz, &fieldp[k][2]);
   }
   // */
} // generated by ComplexKernelBuilder (ck.py) 1.5.1


void ufield_ewald_real_cu(const real (*uind)[3], const real (*uinp)[3],
                          real (*field)[3], real (*fieldp)[3])
{
   const auto& st = *mspatial_v2_unit;
   const real off = switch_off(switch_ewald);


   const PMEUnit pu = ppme_unit;
   const real aewald = pu->aewald;


   int ngrid = get_grid_size(BLOCK_DIM);
   ufield_cu1<EWALD><<<ngrid, BLOCK_DIM, 0, nonblk>>>(
      st.n, TINKER_IMAGE_ARGS, off, st.si1.bit0, st.si2.bit0, st.si3.bit0,
      nmdpuexclude, mdpuexclude, mdpuexclude_scale, st.x, st.y, st.z, st.sorted,
      st.nakpl, st.iakpl, st.niak, st.iak, st.lst, uind, uinp, field, fieldp,
      thole, pdamp, aewald);
}


void ufield_nonewald_cu(const real (*uind)[3], const real (*uinp)[3],
                        real (*field)[3], real (*fieldp)[3])
{
   const auto& st = *mspatial_v2_unit;
   const real off = switch_off(switch_mpole);


   darray::zero(PROCEED_NEW_Q, n, field, fieldp);
   int ngrid = get_grid_size(BLOCK_DIM);
   ufield_cu1<NON_EWALD><<<ngrid, BLOCK_DIM, 0, nonblk>>>(
      st.n, TINKER_IMAGE_ARGS, off, st.si1.bit0, st.si2.bit0, st.si3.bit0,
      nmdpuexclude, mdpuexclude, mdpuexclude_scale, st.x, st.y, st.z, st.sorted,
      st.nakpl, st.iakpl, st.niak, st.iak, st.lst, uind, uinp, field, fieldp,
      thole, pdamp, 0);
}
}
