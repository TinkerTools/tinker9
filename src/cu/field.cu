#include "add.h"
#include "ff/amoeba/elecamoeba.h"
#include "ff/amoeba/epolar.h"
#include "ff/energy.h"
#include "ff/image.h"
#include "ff/pme.h"
#include "ff/spatial.h"
#include "ff/switch.h"
#include "launch.h"
#include "seq/pair_field.h"
#include "seq/triangle.h"
#include "tool/cudalib.h"
#include "tool/gpucard.h"

namespace tinker {
// ck.py Version 2.0.2
template <class ETYP>
__global__
void dfield_cu1(int n, TINKER_IMAGE_PARAMS, real off, const unsigned* restrict dpinfo, int nexclude,
   const int (*restrict exclude)[2], const real (*restrict exclude_scale)[2],
   const real* restrict x, const real* restrict y, const real* restrict z,
   const Spatial::SortedAtom* restrict sorted, int nakpl, const int* restrict iakpl, int niak,
   const int* restrict iak, const int* restrict lst, real (*restrict field)[3],
   real (*restrict fieldp)[3], const real (*restrict rpole)[10], const real* restrict thole,
   const real* restrict pdamp, real aewald)
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
         pair_dfield_v2<ETYP>(r2, xr, yr, zr, scalea, scaleb, aewald, ci[klane], dix[klane],
            diy[klane], diz[klane], qixx[klane], qixy[klane], qixz[klane], qiyy[klane], qiyz[klane],
            qizz[klane], pdi[klane], pti[klane], ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz,
            qkzz, pdk, ptk, fidx, fidy, fidz, fipx, fipy, fipz, fkdx, fkdy, fkdz, fkpx, fkpy, fkpz);
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
            pair_dfield_v2<ETYP>(r2, xr, yr, zr, scalea, scaleb, aewald, ci[klane], dix[klane],
               diy[klane], diz[klane], qixx[klane], qixy[klane], qixz[klane], qiyy[klane],
               qiyz[klane], qizz[klane], pdi[klane], pti[klane], ck, dkx, dky, dkz, qkxx, qkxy,
               qkxz, qkyy, qkyz, qkzz, pdk, ptk, fidx, fidy, fidz, fipx, fipy, fipz, fkdx, fkdy,
               fkdz, fkpx, fkpy, fkpz);
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
            pair_dfield_v2<ETYP>(r2, xr, yr, zr, scalea, scaleb, aewald, ci[klane], dix[klane],
               diy[klane], diz[klane], qixx[klane], qixy[klane], qixz[klane], qiyy[klane],
               qiyz[klane], qizz[klane], pdi[klane], pti[klane], ck, dkx, dky, dkz, qkxx, qkxy,
               qkxz, qkyy, qkyz, qkzz, pdk, ptk, fidx, fidy, fidz, fipx, fipy, fipz, fkdx, fkdy,
               fkdz, fkpx, fkpy, fkpz);
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
   const real off = switchOff(Switch::EWALD);

   const PMEUnit pu = ppme_unit;
   const real aewald = pu->aewald;

   int ngrid = gpuGridSize(BLOCK_DIM);
   ngrid *= BLOCK_DIM;
   int nparallel = std::max(st.niak, st.nakpl) * WARP_SIZE;
   nparallel = std::max(nparallel, ngrid);
   launch_k1s(g::s0, nparallel, dfield_cu1<EWALD>, //
      st.n, TINKER_IMAGE_ARGS, off, st.si3.bit0, ndpexclude, dpexclude, dpexclude_scale, st.x, st.y,
      st.z, st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst, field, fieldp, rpole, thole,
      pdamp, aewald);
}

void dfield_nonewald_cu(real (*field)[3], real (*fieldp)[3])
{
   const auto& st = *mspatial_v2_unit;
   const real off = switchOff(Switch::MPOLE);

   darray::zero(g::q0, n, field, fieldp);
   int ngrid = gpuGridSize(BLOCK_DIM);
   ngrid *= BLOCK_DIM;
   int nparallel = std::max(st.niak, st.nakpl) * WARP_SIZE;
   nparallel = std::max(nparallel, ngrid);
   launch_k1s(g::s0, nparallel, dfield_cu1<NON_EWALD>, //
      st.n, TINKER_IMAGE_ARGS, off, st.si3.bit0, ndpexclude, dpexclude, dpexclude_scale, st.x, st.y,
      st.z, st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst, field, fieldp, rpole, thole,
      pdamp, 0);
}

// ck.py Version 2.0.2
template <class ETYP>
__global__
void ufield_cu1(int n, TINKER_IMAGE_PARAMS, real off, const unsigned* restrict uinfo, int nexclude,
   const int (*restrict exclude)[2], const real* restrict exclude_scale, const real* restrict x,
   const real* restrict y, const real* restrict z, const Spatial::SortedAtom* restrict sorted,
   int nakpl, const int* restrict iakpl, int niak, const int* restrict iak, const int* restrict lst,
   const real (*restrict uind)[3], const real (*restrict uinp)[3], real (*restrict field)[3],
   real (*restrict fieldp)[3], const real* restrict thole, const real* restrict pdamp, real aewald)
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
   __shared__ real uidx[BLOCK_DIM];
   __shared__ real uidy[BLOCK_DIM];
   __shared__ real uidz[BLOCK_DIM];
   __shared__ real uipx[BLOCK_DIM];
   __shared__ real uipy[BLOCK_DIM];
   __shared__ real uipz[BLOCK_DIM];
   __shared__ real pdi[BLOCK_DIM];
   __shared__ real pti[BLOCK_DIM];
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
      real scalea = exclude_scale[ii];

      xi = x[i];
      yi = y[i];
      zi = z[i];
      xk = x[k];
      yk = y[k];
      zk = z[k];
      uidx[klane] = uind[i][0];
      uidy[klane] = uind[i][1];
      uidz[klane] = uind[i][2];
      uipx[klane] = uinp[i][0];
      uipy[klane] = uinp[i][1];
      uipz[klane] = uinp[i][2];
      pdi[klane] = pdamp[i];
      pti[klane] = thole[i];
      ukdx = uind[k][0];
      ukdy = uind[k][1];
      ukdz = uind[k][2];
      ukpx = uinp[k][0];
      ukpy = uinp[k][1];
      ukpz = uinp[k][2];
      pdk = pdamp[k];
      ptk = thole[k];

      constexpr bool incl = true;
      real xr = xk - xi;
      real yr = yk - yi;
      real zr = zk - zi;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         pair_ufield_v2<ETYP>(r2, xr, yr, zr, scalea, aewald, uidx[klane], uidy[klane], uidz[klane],
            uipx[klane], uipy[klane], uipz[klane], pdi[klane], pti[klane], ukdx, ukdy, ukdz, ukpx,
            ukpy, ukpz, pdk, ptk, fidx, fidy, fidz, fipx, fipy, fipz, fkdx, fkdy, fkdz, fkpx, fkpy,
            fkpz);
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

      uidx[threadIdx.x] = uind[i][0];
      uidy[threadIdx.x] = uind[i][1];
      uidz[threadIdx.x] = uind[i][2];
      uipx[threadIdx.x] = uinp[i][0];
      uipy[threadIdx.x] = uinp[i][1];
      uipz[threadIdx.x] = uinp[i][2];
      pdi[threadIdx.x] = pdamp[i];
      pti[threadIdx.x] = thole[i];
      ukdx = uind[k][0];
      ukdy = uind[k][1];
      ukdz = uind[k][2];
      ukpx = uinp[k][0];
      ukpy = uinp[k][1];
      ukpz = uinp[k][2];
      pdk = pdamp[k];
      ptk = thole[k];

      unsigned int uinfo0 = uinfo[iw * WARP_SIZE + ilane];
      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         bool incl = iid < kid and kid < n;
         int srcmask = 1 << srclane;
         incl = incl and (uinfo0 & srcmask) == 0;
         real scalea = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_ufield_v2<ETYP>(r2, xr, yr, zr, scalea, aewald, uidx[klane], uidy[klane],
               uidz[klane], uipx[klane], uipy[klane], uipz[klane], pdi[klane], pti[klane], ukdx,
               ukdy, ukdz, ukpx, ukpy, ukpz, pdk, ptk, fidx, fidy, fidz, fipx, fipy, fipz, fkdx,
               fkdy, fkdz, fkpx, fkpy, fkpz);
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

      uidx[threadIdx.x] = uind[i][0];
      uidy[threadIdx.x] = uind[i][1];
      uidz[threadIdx.x] = uind[i][2];
      uipx[threadIdx.x] = uinp[i][0];
      uipy[threadIdx.x] = uinp[i][1];
      uipz[threadIdx.x] = uinp[i][2];
      pdi[threadIdx.x] = pdamp[i];
      pti[threadIdx.x] = thole[i];
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
         real scalea = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_ufield_v2<ETYP>(r2, xr, yr, zr, scalea, aewald, uidx[klane], uidy[klane],
               uidz[klane], uipx[klane], uipy[klane], uipz[klane], pdi[klane], pti[klane], ukdx,
               ukdy, ukdz, ukpx, ukpy, ukpz, pdk, ptk, fidx, fidy, fidz, fipx, fipy, fipz, fkdx,
               fkdy, fkdz, fkpx, fkpy, fkpz);
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

void ufield_ewald_real_cu(
   const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3])
{
   const auto& st = *mspatial_v2_unit;
   const real off = switchOff(Switch::EWALD);

   const PMEUnit pu = ppme_unit;
   const real aewald = pu->aewald;

   int ngrid = gpuGridSize(BLOCK_DIM);
   ngrid *= BLOCK_DIM;
   int nparallel = std::max(st.niak, st.nakpl) * WARP_SIZE;
   nparallel = std::max(nparallel, ngrid);
   launch_k1s(g::s0, nparallel, ufield_cu1<EWALD>, //
      st.n, TINKER_IMAGE_ARGS, off, st.si4.bit0, nuexclude, uexclude, uexclude_scale, st.x, st.y,
      st.z, st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst, uind, uinp, field, fieldp,
      thole, pdamp, aewald);
}

void ufield_nonewald_cu(
   const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3])
{
   const auto& st = *mspatial_v2_unit;
   const real off = switchOff(Switch::MPOLE);

   darray::zero(g::q0, n, field, fieldp);
   int ngrid = gpuGridSize(BLOCK_DIM);
   ngrid *= BLOCK_DIM;
   int nparallel = std::max(st.niak, st.nakpl) * WARP_SIZE;
   nparallel = std::max(nparallel, ngrid);
   launch_k1s(g::s0, nparallel, ufield_cu1<NON_EWALD>, //
      st.n, TINKER_IMAGE_ARGS, off, st.si4.bit0, nuexclude, uexclude, uexclude_scale, st.x, st.y,
      st.z, st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst, uind, uinp, field, fieldp,
      thole, pdamp, 0);
}
}
