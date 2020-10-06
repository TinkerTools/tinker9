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
template <class ETYP>
__launch_bounds__(BLOCK_DIM) __global__
void dfield_cu1(int n, TINKER_IMAGE_PARAMS, real off,
                const unsigned* restrict minfo, const unsigned* restrict dpinfo,
                const unsigned* restrict uinfo, int nexclude,
                const int (*restrict exclude)[2],
                const real (*restrict exclude_scale)[4], const real* restrict x,
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
      pdk = pdamp[k];
      ptk = thole[k];


      int klane = threadIdx.x;
      constexpr bool incl = true;
      real xr = xk - xi;
      real yr = yk - yi;
      real zr = zk - zi;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         pair_dfield_v2<ETYP>(
            r2, xr, yr, zr, scaleb, scalec, aewald, ci, dix, diy, diz, qixx,
            qixy, qixz, qiyy, qiyz, qizz, pdi, pti, ck, dkx, dky, dkz, qkxx,
            qkxy, qkxz, qkyy, qkyz, qkzz, pdk, ptk, shfidx[klane],
            shfidy[klane], shfidz[klane], shfipx[klane], shfipy[klane],
            shfipz[klane], fkdx, fkdy, fkdz, fkpx, fkpy, fkpz);
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
         real pdi = shpdi[klane];
         real pti = shpti[klane];


         bool incl = iid < kid and kid < n;
         incl = incl and (minfo0 & srcmask) == 0 and
            (dpinfo0 & srcmask) == 0 and (uinfo0 & srcmask) == 0;
         real scaleb = 1;
         real scalec = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_dfield_v2<ETYP>(
               r2, xr, yr, zr, scaleb, scalec, aewald, ci, dix, diy, diz, qixx,
               qixy, qixz, qiyy, qiyz, qizz, pdi, pti, ck, dkx, dky, dkz, qkxx,
               qkxy, qkxz, qkyy, qkyz, qkzz, pdk, ptk, shfidx[klane],
               shfidy[klane], shfidz[klane], shfipx[klane], shfipy[klane],
               shfipz[klane], fkdx, fkdy, fkdz, fkpx, fkpy, fkpz);
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
         real pdi = shpdi[klane];
         real pti = shpti[klane];


         bool incl = atomk > 0;
         real scaleb = 1;
         real scalec = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_dfield_v2<ETYP>(
               r2, xr, yr, zr, scaleb, scalec, aewald, ci, dix, diy, diz, qixx,
               qixy, qixz, qiyy, qiyz, qizz, pdi, pti, ck, dkx, dky, dkz, qkxx,
               qkxy, qkxz, qkyy, qkyz, qkzz, pdk, ptk, shfidx[klane],
               shfidy[klane], shfidz[klane], shfipx[klane], shfipy[klane],
               shfipz[klane], fkdx, fkdy, fkdz, fkpx, fkpy, fkpz);
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


void dfield_ewald_real_cu(real (*field)[3], real (*fieldp)[3])
{
   const auto& st = *mspatial_v2_unit;
   const real off = switch_off(switch_ewald);


   const PMEUnit pu = ppme_unit;
   const real aewald = pu->aewald;


   int ngrid = get_grid_size(BLOCK_DIM);
   dfield_cu1<EWALD><<<ngrid, BLOCK_DIM, 0, nonblk>>>(
      st.n, TINKER_IMAGE_ARGS, off, st.si1.bit0, st.si2.bit0, st.si3.bit0,
      nmdpuexclude, mdpuexclude, mdpuexclude_scale, st.x, st.y, st.z, st.sorted,
      st.nakpl, st.iakpl, st.niak, st.iak, st.lst, field, fieldp, rpole, thole,
      pdamp, aewald);
}


void dfield_nonewald_cu(real (*field)[3], real (*fieldp)[3])
{
   const auto& st = *mspatial_v2_unit;
   const real off = switch_off(switch_mpole);


   darray::zero(PROCEED_NEW_Q, n, field, fieldp);
   int ngrid = get_grid_size(BLOCK_DIM);
   dfield_cu1<NON_EWALD><<<ngrid, BLOCK_DIM, 0, nonblk>>>(
      st.n, TINKER_IMAGE_ARGS, off, st.si1.bit0, st.si2.bit0, st.si3.bit0,
      nmdpuexclude, mdpuexclude, mdpuexclude_scale, st.x, st.y, st.z, st.sorted,
      st.nakpl, st.iakpl, st.niak, st.iak, st.lst, field, fieldp, rpole, thole,
      pdamp, 0);
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
