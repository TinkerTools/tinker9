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
// ck.py Version 2.0.2
template <class ETYP>
__global__
void dfield_cu1(int n, TINKER_IMAGE_PARAMS, real off,
                const unsigned* restrict dpinfo, int nexclude,
                const int (*restrict exclude)[2],
                const real (*restrict exclude_scale)[2], const real* restrict x,
                const real* restrict y, const real* restrict z,
                const Spatial::SortedAtom* restrict sorted, int nakpl,
                const int* restrict iakpl, int niak, const int* restrict iak,
                const int* restrict lst, real (*restrict field)[3],
                real (*restrict fieldp)[3], const real (*restrict rpole)[10],
                const real* restrict thole, const real* restrict pdamp,
                real aewald)
{
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);


   real xi;
   real yi;
   real zi;
   real xk;
   real yk;
   real zk;
   real fidx;
   real fidy;
   real fidz;
   real fipx;
   real fipy;
   real fipz;
   real fkdx;
   real fkdy;
   real fkdz;
   real fkpx;
   real fkpy;
   real fkpz;
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
   real pdk;
   real ptk;


   //* /
   for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
      const int klane = threadIdx.x;
      fidx = 0;
      fidy = 0;
      fidz = 0;
      fipx = 0;
      fipy = 0;
      fipz = 0;
      fkdx = 0;
      fkdy = 0;
      fkdz = 0;
      fkpx = 0;
      fkpy = 0;
      fkpz = 0;


      int i = exclude[ii][0];
      int k = exclude[ii][1];
      real scalea = exclude_scale[ii][0];
      real scaleb = exclude_scale[ii][1];


      xi = x[i];
      yi = y[i];
      zi = z[i];
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
      pdk = pdamp[k];
      ptk = thole[k];


      constexpr bool incl = true;
      real xr = xk - xi;
      real yr = yk - yi;
      real zr = zk - zi;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         pair_dfield_v2<ETYP>(
            r2, xr, yr, zr, scalea, scaleb, aewald, ci[klane], dix[klane],
            diy[klane], diz[klane], qixx[klane], qixy[klane], qixz[klane],
            qiyy[klane], qiyz[klane], qizz[klane], pdi[klane], pti[klane], ck,
            dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, pdk, ptk, fidx,
            fidy, fidz, fipx, fipy, fipz, fkdx, fkdy, fkdz, fkpx, fkpy, fkpz);
      } // end if (include)


      atomic_add(fidx, &field[i][0]);
      atomic_add(fidy, &field[i][1]);
      atomic_add(fidz, &field[i][2]);
      atomic_add(fipx, &fieldp[i][0]);
      atomic_add(fipy, &fieldp[i][1]);
      atomic_add(fipz, &fieldp[i][2]);
      atomic_add(fkdx, &field[k][0]);
      atomic_add(fkdy, &field[k][1]);
      atomic_add(fkdz, &field[k][2]);
      atomic_add(fkpx, &fieldp[k][0]);
      atomic_add(fkpy, &fieldp[k][1]);
      atomic_add(fkpz, &fieldp[k][2]);
   }
   // */


   for (int iw = iwarp; iw < nakpl; iw += nwarp) {
      fidx = 0;
      fidy = 0;
      fidz = 0;
      fipx = 0;
      fipy = 0;
      fipz = 0;
      fkdx = 0;
      fkdy = 0;
      fkdz = 0;
      fkpx = 0;
      fkpy = 0;
      fkpz = 0;


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
      pdk = pdamp[k];
      ptk = thole[k];


      unsigned int dpinfo0 = dpinfo[iw * WARP_SIZE + ilane];
      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         bool incl = iid < kid and kid < n;
         int srcmask = 1 << srclane;
         incl = incl and (dpinfo0 & srcmask) == 0;
         real scalea = 1;
         real scaleb = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_dfield_v2<ETYP>(
               r2, xr, yr, zr, scalea, scaleb, aewald, ci[klane], dix[klane],
               diy[klane], diz[klane], qixx[klane], qixy[klane], qixz[klane],
               qiyy[klane], qiyz[klane], qizz[klane], pdi[klane], pti[klane],
               ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, pdk, ptk,
               fidx, fidy, fidz, fipx, fipy, fipz, fkdx, fkdy, fkdz, fkpx, fkpy,
               fkpz);
         } // end if (include)


         iid = __shfl_sync(ALL_LANES, iid, ilane + 1);
         xi = __shfl_sync(ALL_LANES, xi, ilane + 1);
         yi = __shfl_sync(ALL_LANES, yi, ilane + 1);
         zi = __shfl_sync(ALL_LANES, zi, ilane + 1);
         fidx = __shfl_sync(ALL_LANES, fidx, ilane + 1);
         fidy = __shfl_sync(ALL_LANES, fidy, ilane + 1);
         fidz = __shfl_sync(ALL_LANES, fidz, ilane + 1);
         fipx = __shfl_sync(ALL_LANES, fipx, ilane + 1);
         fipy = __shfl_sync(ALL_LANES, fipy, ilane + 1);
         fipz = __shfl_sync(ALL_LANES, fipz, ilane + 1);
      }


      atomic_add(fidx, &field[i][0]);
      atomic_add(fidy, &field[i][1]);
      atomic_add(fidz, &field[i][2]);
      atomic_add(fipx, &fieldp[i][0]);
      atomic_add(fipy, &fieldp[i][1]);
      atomic_add(fipz, &fieldp[i][2]);
      atomic_add(fkdx, &field[k][0]);
      atomic_add(fkdy, &field[k][1]);
      atomic_add(fkdz, &field[k][2]);
      atomic_add(fkpx, &fieldp[k][0]);
      atomic_add(fkpy, &fieldp[k][1]);
      atomic_add(fkpz, &fieldp[k][2]);
   }


   for (int iw = iwarp; iw < niak; iw += nwarp) {
      fidx = 0;
      fidy = 0;
      fidz = 0;
      fipx = 0;
      fipy = 0;
      fipz = 0;
      fkdx = 0;
      fkdy = 0;
      fkdz = 0;
      fkpx = 0;
      fkpy = 0;
      fkpz = 0;


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
      pdk = pdamp[k];
      ptk = thole[k];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         bool incl = atomk > 0;
         real scalea = 1;
         real scaleb = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_dfield_v2<ETYP>(
               r2, xr, yr, zr, scalea, scaleb, aewald, ci[klane], dix[klane],
               diy[klane], diz[klane], qixx[klane], qixy[klane], qixz[klane],
               qiyy[klane], qiyz[klane], qizz[klane], pdi[klane], pti[klane],
               ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, pdk, ptk,
               fidx, fidy, fidz, fipx, fipy, fipz, fkdx, fkdy, fkdz, fkpx, fkpy,
               fkpz);
         } // end if (include)


         xi = __shfl_sync(ALL_LANES, xi, ilane + 1);
         yi = __shfl_sync(ALL_LANES, yi, ilane + 1);
         zi = __shfl_sync(ALL_LANES, zi, ilane + 1);
         fidx = __shfl_sync(ALL_LANES, fidx, ilane + 1);
         fidy = __shfl_sync(ALL_LANES, fidy, ilane + 1);
         fidz = __shfl_sync(ALL_LANES, fidz, ilane + 1);
         fipx = __shfl_sync(ALL_LANES, fipx, ilane + 1);
         fipy = __shfl_sync(ALL_LANES, fipy, ilane + 1);
         fipz = __shfl_sync(ALL_LANES, fipz, ilane + 1);
      }


      atomic_add(fidx, &field[i][0]);
      atomic_add(fidy, &field[i][1]);
      atomic_add(fidz, &field[i][2]);
      atomic_add(fipx, &fieldp[i][0]);
      atomic_add(fipy, &fieldp[i][1]);
      atomic_add(fipz, &fieldp[i][2]);
      atomic_add(fkdx, &field[k][0]);
      atomic_add(fkdy, &field[k][1]);
      atomic_add(fkdz, &field[k][2]);
      atomic_add(fkpx, &fieldp[k][0]);
      atomic_add(fkpy, &fieldp[k][1]);
      atomic_add(fkpz, &fieldp[k][2]);
   }
}


void dfield_ewald_real_cu(real (*field)[3], real (*fieldp)[3])
{
   const auto& st = *mspatial_v2_unit;
   const real off = switch_off(switch_ewald);


   const PMEUnit pu = ppme_unit;
   const real aewald = pu->aewald;


   int ngrid = get_grid_size(BLOCK_DIM);
   dfield_cu1<EWALD><<<ngrid, BLOCK_DIM, 0, nonblk>>>(
      st.n, TINKER_IMAGE_ARGS, off, st.si3.bit0, ndpexclude, dpexclude,
      dpexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl, st.niak,
      st.iak, st.lst, field, fieldp, rpole, thole, pdamp, aewald);
}


void dfield_nonewald_cu(real (*field)[3], real (*fieldp)[3])
{
   const auto& st = *mspatial_v2_unit;
   const real off = switch_off(switch_mpole);


   darray::zero(PROCEED_NEW_Q, n, field, fieldp);
   int ngrid = get_grid_size(BLOCK_DIM);
   dfield_cu1<NON_EWALD><<<ngrid, BLOCK_DIM, 0, nonblk>>>(
      st.n, TINKER_IMAGE_ARGS, off, st.si3.bit0, ndpexclude, dpexclude,
      dpexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl, st.niak,
      st.iak, st.lst, field, fieldp, rpole, thole, pdamp, 0);
}


template <class ETYP>
__launch_bounds__(BLOCK_DIM) __global__
void ufield_cu1(int n, TINKER_IMAGE_PARAMS, real off,
                const unsigned* restrict mdpuinfo, int nexclude,
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


      unsigned int mdpuinfo0 = mdpuinfo[iw * WARP_SIZE + ilane];


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
         incl = incl and (mdpuinfo0 & srcmask) == 0;
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
      st.n, TINKER_IMAGE_ARGS, off, st.si1.bit0, nmdpuexclude, mdpuexclude,
      mdpuexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl,
      st.niak, st.iak, st.lst, uind, uinp, field, fieldp, thole, pdamp, aewald);
}


void ufield_nonewald_cu(const real (*uind)[3], const real (*uinp)[3],
                        real (*field)[3], real (*fieldp)[3])
{
   const auto& st = *mspatial_v2_unit;
   const real off = switch_off(switch_mpole);


   darray::zero(PROCEED_NEW_Q, n, field, fieldp);
   int ngrid = get_grid_size(BLOCK_DIM);
   ufield_cu1<NON_EWALD><<<ngrid, BLOCK_DIM, 0, nonblk>>>(
      st.n, TINKER_IMAGE_ARGS, off, st.si1.bit0, nmdpuexclude, mdpuexclude,
      mdpuexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl,
      st.niak, st.iak, st.lst, uind, uinp, field, fieldp, thole, pdamp, 0);
}
}
