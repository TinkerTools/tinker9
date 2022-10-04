#include "ff/amoebamod.h"
#include "ff/image.h"
#include "ff/pme.h"
#include "ff/spatial.h"
#include "ff/switch.h"
#include "seq/launch.h"
#include "seq/pair_field.h"
#include "seq/triangle.h"

#define TINKER9_POLPAIR 2
#include "ff/amoebacumod.h"

namespace tinker {
__global__
void dfieldEwaldRecipSelfP2_cu1(int n, real (*restrict field)[3], real term,
   const real (*restrict rpole)[MPL_TOTAL], const real (*restrict cphi)[10])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      real dix = rpole[i][MPL_PME_X];
      real diy = rpole[i][MPL_PME_Y];
      real diz = rpole[i][MPL_PME_Z];
      field[i][0] += (-cphi[i][1] + term * dix);
      field[i][1] += (-cphi[i][2] + term * diy);
      field[i][2] += (-cphi[i][3] + term * diz);
   }
}

void dfieldEwaldRecipSelfP2_cu(real (*field)[3])
{
   const PMEUnit pu = ppme_unit;
   const real aewald = pu->aewald;
   const real term = aewald * aewald * aewald * 4 / 3 / sqrtpi;

   launch_k1s(g::s0, n, dfieldEwaldRecipSelfP2_cu1, n, field, term, rpole,
      cphi);
}

#if TINKER9_POLPAIR == 0
// ck.py Version 2.0.2
template <class ETYP>
__global__
void dfield_cu1(int n, TINKER_IMAGE_PARAMS, real off,
   const unsigned* restrict dpinfo, int nexclude,
   const int (*restrict exclude)[2], const real (*restrict exclude_scale)[2],
   const real* restrict x, const real* restrict y, const real* restrict z,
   const Spatial::SortedAtom* restrict sorted, int nakpl,
   const int* restrict iakpl, int niak, const int* restrict iak,
   const int* restrict lst, real (*restrict field)[3],
   real (*restrict fieldp)[3],
#   if TINKER9_POLPAIR == 0
   const real (*restrict rpole)[10], const real* restrict thole,
   const real* restrict pdamp,
#   endif
   real aewald)
{
#   if TINKER9_POLPAIR
   using d::pdamp;
   using d::rpole;
#      if TINKER9_POLPAIR == 1
   using d::thole;
#      elif TINKER9_POLPAIR == 2
   using d::jpolar;
   using d::njpolar;
   using d::thlval;
#      endif
#   endif

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
#   if TINKER9_POLPAIR == 1
   __shared__ real pti[BLOCK_DIM];
#   elif TINKER9_POLPAIR == 2
   __shared__ int jpi[BLOCK_DIM];
#   endif
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
#   if TINKER9_POLPAIR == 1
   real ptk;
#   elif TINKER9_POLPAIR == 2
   int jpk;
#   endif

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
      pdk = pdamp[k];
      ptk = thole[k];

      constexpr bool incl = true;
      real xr = xk - xi;
      real yr = yk - yi;
      real zr = zk - zi;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         pair_dfield_v2<ETYP>(r2, xr, yr, zr, scalea, scaleb, aewald, ci[klane],
            dix[klane], diy[klane], diz[klane], qixx[klane], qixy[klane],
            qixz[klane], qiyy[klane], qiyz[klane], qizz[klane], pdi[klane],
            pti[klane], ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz,
            pdk, ptk, fidx, fidy, fidz, fipx, fipy, fipz, fkdx, fkdy, fkdz,
            fkpx, fkpy, fkpz);
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
            pair_dfield_v2<ETYP>(r2, xr, yr, zr, scalea, scaleb, aewald,
               ci[klane], dix[klane], diy[klane], diz[klane], qixx[klane],
               qixy[klane], qixz[klane], qiyy[klane], qiyz[klane], qizz[klane],
               pdi[klane], pti[klane], ck, dkx, dky, dkz, qkxx, qkxy, qkxz,
               qkyy, qkyz, qkzz, pdk, ptk, fidx, fidy, fidz, fipx, fipy, fipz,
               fkdx, fkdy, fkdz, fkpx, fkpy, fkpz);
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
            pair_dfield_v2<ETYP>(r2, xr, yr, zr, scalea, scaleb, aewald,
               ci[klane], dix[klane], diy[klane], diz[klane], qixx[klane],
               qixy[klane], qixz[klane], qiyy[klane], qiyz[klane], qizz[klane],
               pdi[klane], pti[klane], ck, dkx, dky, dkz, qkxx, qkxy, qkxz,
               qkyy, qkyz, qkzz, pdk, ptk, fidx, fidy, fidz, fipx, fipy, fipz,
               fkdx, fkdy, fkdz, fkpx, fkpy, fkpz);
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

void dfieldEwaldReal_cu(real (*field)[3], real (*fieldp)[3])
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
      st.n, TINKER_IMAGE_ARGS, off, st.si3.bit0, ndpexclude, dpexclude,
      dpexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl, st.niak,
      st.iak, st.lst, field, fieldp,
#   if TINKER9_POLPAIR
#   else
      rpole, thole, pdamp,
#   endif
      aewald);
}

void dfieldNonEwald_cu(real (*field)[3], real (*fieldp)[3])
{
   const auto& st = *mspatial_v2_unit;
   const real off = switchOff(Switch::MPOLE);

   darray::zero(g::q0, n, field, fieldp);
   int ngrid = gpuGridSize(BLOCK_DIM);
   ngrid *= BLOCK_DIM;
   int nparallel = std::max(st.niak, st.nakpl) * WARP_SIZE;
   nparallel = std::max(nparallel, ngrid);
   launch_k1s(g::s0, nparallel, dfield_cu1<NON_EWALD>, //
      st.n, TINKER_IMAGE_ARGS, off, st.si3.bit0, ndpexclude, dpexclude,
      dpexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl, st.niak,
      st.iak, st.lst, field, fieldp,
#   if TINKER9_POLPAIR
#   else
      rpole, thole, pdamp,
#   endif
      0);
}
#elif TINKER9_POLPAIR == 2
// ck.py Version 2.0.2
template <class ETYP>
__global__
void dfield_cu1(int n, TINKER_IMAGE_PARAMS, real off,
   const unsigned* restrict dpinfo, int nexclude,
   const int (*restrict exclude)[2], const real (*restrict exclude_scale)[2],
   const real* restrict x, const real* restrict y, const real* restrict z,
   const Spatial::SortedAtom* restrict sorted, int nakpl,
   const int* restrict iakpl, int niak, const int* restrict iak,
   const int* restrict lst, real (*restrict field)[3],
   real (*restrict fieldp)[3], real aewald)
{

   using d::jpolar;
   using d::njpolar;
   using d::pdamp;
   using d::rpole;
   using d::thlval;

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
   __shared__ int jpi[BLOCK_DIM];
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
   int jpk;

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
      pdi[klane] = pdamp[i];
      jpi[klane] = jpolar[i];
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
      pdk = pdamp[k];
      jpk = jpolar[k];

      constexpr bool incl = true;
      real xr = xk - xi;
      real yr = yk - yi;
      real zr = zk - zi;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         real pga = thlval[njpolar * jpi[klane] + jpk];
         pair_dfield_v2<ETYP>(r2, xr, yr, zr, scalea, scaleb, aewald, ci[klane],
            dix[klane], diy[klane], diz[klane], qixx[klane], qixy[klane],
            qixz[klane], qiyy[klane], qiyz[klane], qizz[klane], pdi[klane], pga,
            ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, pdk, pga,
            fidx, fidy, fidz, fipx, fipy, fipz, fkdx, fkdy, fkdz, fkpx, fkpy,
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
      pdi[threadIdx.x] = pdamp[i];
      jpi[threadIdx.x] = jpolar[i];
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
      pdk = pdamp[k];
      jpk = jpolar[k];
      __syncwarp();

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
            real pga = thlval[njpolar * jpi[klane] + jpk];
            pair_dfield_v2<ETYP>(r2, xr, yr, zr, scalea, scaleb, aewald,
               ci[klane], dix[klane], diy[klane], diz[klane], qixx[klane],
               qixy[klane], qixz[klane], qiyy[klane], qiyz[klane], qizz[klane],
               pdi[klane], pga, ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz,
               qkzz, pdk, pga, fidx, fidy, fidz, fipx, fipy, fipz, fkdx, fkdy,
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
   __syncwarp();

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
      pdi[threadIdx.x] = pdamp[i];
      jpi[threadIdx.x] = jpolar[i];
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
      pdk = pdamp[k];
      jpk = jpolar[k];
      __syncwarp();

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
            real pga = thlval[njpolar * jpi[klane] + jpk];
            pair_dfield_v2<ETYP>(r2, xr, yr, zr, scalea, scaleb, aewald,
               ci[klane], dix[klane], diy[klane], diz[klane], qixx[klane],
               qixy[klane], qixz[klane], qiyy[klane], qiyz[klane], qizz[klane],
               pdi[klane], pga, ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz,
               qkzz, pdk, pga, fidx, fidy, fidz, fipx, fipy, fipz, fkdx, fkdy,
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
   __syncwarp();
}

void dfieldEwaldReal_cu(real (*field)[3], real (*fieldp)[3])
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
      st.n, TINKER_IMAGE_ARGS, off, st.si3.bit0, ndpexclude, dpexclude,
      dpexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl, st.niak,
      st.iak, st.lst, field, fieldp, aewald);
}

void dfieldNonEwald_cu(real (*field)[3], real (*fieldp)[3])
{
   const auto& st = *mspatial_v2_unit;
   const real off = switchOff(Switch::MPOLE);

   darray::zero(g::q0, n, field, fieldp);
   int ngrid = gpuGridSize(BLOCK_DIM);
   ngrid *= BLOCK_DIM;
   int nparallel = std::max(st.niak, st.nakpl) * WARP_SIZE;
   nparallel = std::max(nparallel, ngrid);
   launch_k1s(g::s0, nparallel, dfield_cu1<NON_EWALD>, //
      st.n, TINKER_IMAGE_ARGS, off, st.si3.bit0, ndpexclude, dpexclude,
      dpexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl, st.niak,
      st.iak, st.lst, field, fieldp, 0);
}
#endif
}

namespace tinker {
__global__
void ufieldEwaldRecipSelfP1_cu1(int n, const real (*restrict uind)[3],
   const real (*restrict uinp)[3], real (*restrict field)[3],
   real (*restrict fieldp)[3], const real (*restrict fdip_phi1)[10],
   const real (*restrict fdip_phi2)[10], real term, int nfft1, int nfft2,
   int nfft3, TINKER_IMAGE_PARAMS)
{
   real a[3][3];
   a[0][0] = nfft1 * recipa.x;
   a[1][0] = nfft2 * recipb.x;
   a[2][0] = nfft3 * recipc.x;
   a[0][1] = nfft1 * recipa.y;
   a[1][1] = nfft2 * recipb.y;
   a[2][1] = nfft3 * recipc.y;
   a[0][2] = nfft1 * recipa.z;
   a[1][2] = nfft2 * recipb.z;
   a[2][2] = nfft3 * recipc.z;

   if (uinp) {
      for (int i = ITHREAD; i < n; i += STRIDE) {
         for (int j = 0; j < 3; ++j) {
            real df1 = a[0][j] * fdip_phi1[i][1] + a[1][j] * fdip_phi1[i][2]
               + a[2][j] * fdip_phi1[i][3];
            real df2 = a[0][j] * fdip_phi2[i][1] + a[1][j] * fdip_phi2[i][2]
               + a[2][j] * fdip_phi2[i][3];
            field[i][j] += (term * uind[i][j] - df1);
            fieldp[i][j] += (term * uinp[i][j] - df2);
         }
      }
   } else {
      for (int i = ITHREAD; i < n; i += STRIDE) {
         for (int j = 0; j < 3; ++j) {
            real df1 = a[0][j] * fdip_phi1[i][1] + a[1][j] * fdip_phi1[i][2]
               + a[2][j] * fdip_phi1[i][3];
            field[i][j] += (term * uind[i][j] - df1);
         }
      }
   }
}

void ufieldEwaldRecipSelfP1_cu(const real (*uind)[3], const real (*uinp)[3], //
   real (*field)[3], real (*fieldp)[3])
{
   const PMEUnit pu = ppme_unit;
   const real aewald = pu->aewald;
   const real term = aewald * aewald * aewald * 4 / 3 / sqrtpi;
   const int nfft1 = pu->nfft1;
   const int nfft2 = pu->nfft2;
   const int nfft3 = pu->nfft3;

   launch_k1s(g::s0, n, ufieldEwaldRecipSelfP1_cu1, n, uind, uinp, field,
      fieldp, fdip_phi1, fdip_phi2, term, nfft1, nfft2, nfft3,
      TINKER_IMAGE_ARGS);
}

#if TINKER9_POLPAIR == 0
// ck.py Version 2.0.2
template <class ETYP>
__global__
void ufield_cu1(int n, TINKER_IMAGE_PARAMS, real off,
   const unsigned* restrict uinfo, int nexclude,
   const int (*restrict exclude)[2], const real* restrict exclude_scale,
   const real* restrict x, const real* restrict y, const real* restrict z,
   const Spatial::SortedAtom* restrict sorted, int nakpl,
   const int* restrict iakpl, int niak, const int* restrict iak,
   const int* restrict lst, const real (*restrict uind)[3],
   const real (*restrict uinp)[3], real (*restrict field)[3],
   real (*restrict fieldp)[3],
#   if TINKER9_POLPAIR
#   else
   const real* restrict thole, const real* restrict pdamp,
#   endif
   real aewald)
{
#   if TINKER9_POLPAIR
   using d::pdamp;
   using d::thole;
#   else
#   endif

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
         pair_ufield_v2<ETYP>(r2, xr, yr, zr, scalea, aewald, uidx[klane],
            uidy[klane], uidz[klane], uipx[klane], uipy[klane], uipz[klane],
            pdi[klane], pti[klane], ukdx, ukdy, ukdz, ukpx, ukpy, ukpz, pdk,
            ptk, fidx, fidy, fidz, fipx, fipy, fipz, fkdx, fkdy, fkdz, fkpx,
            fkpy, fkpz);
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
            pair_ufield_v2<ETYP>(r2, xr, yr, zr, scalea, aewald, uidx[klane],
               uidy[klane], uidz[klane], uipx[klane], uipy[klane], uipz[klane],
               pdi[klane], pti[klane], ukdx, ukdy, ukdz, ukpx, ukpy, ukpz, pdk,
               ptk, fidx, fidy, fidz, fipx, fipy, fipz, fkdx, fkdy, fkdz, fkpx,
               fkpy, fkpz);
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
            pair_ufield_v2<ETYP>(r2, xr, yr, zr, scalea, aewald, uidx[klane],
               uidy[klane], uidz[klane], uipx[klane], uipy[klane], uipz[klane],
               pdi[klane], pti[klane], ukdx, ukdy, ukdz, ukpx, ukpy, ukpz, pdk,
               ptk, fidx, fidy, fidz, fipx, fipy, fipz, fkdx, fkdy, fkdz, fkpx,
               fkpy, fkpz);
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

void ufieldEwaldReal_cu(const real (*uind)[3], const real (*uinp)[3],
   real (*field)[3], real (*fieldp)[3])
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
      st.n, TINKER_IMAGE_ARGS, off, st.si4.bit0, nuexclude, uexclude,
      uexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl, st.niak,
      st.iak, st.lst, uind, uinp, field, fieldp,
#   if TINKER9_POLPAIR
#   else
      thole, pdamp,
#   endif
      aewald);
}

void ufieldNonEwald_cu(const real (*uind)[3], const real (*uinp)[3],
   real (*field)[3], real (*fieldp)[3])
{
   const auto& st = *mspatial_v2_unit;
   const real off = switchOff(Switch::MPOLE);

   darray::zero(g::q0, n, field, fieldp);
   int ngrid = gpuGridSize(BLOCK_DIM);
   ngrid *= BLOCK_DIM;
   int nparallel = std::max(st.niak, st.nakpl) * WARP_SIZE;
   nparallel = std::max(nparallel, ngrid);
   launch_k1s(g::s0, nparallel, ufield_cu1<NON_EWALD>, //
      st.n, TINKER_IMAGE_ARGS, off, st.si4.bit0, nuexclude, uexclude,
      uexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl, st.niak,
      st.iak, st.lst, uind, uinp, field, fieldp,
#   if TINKER9_POLPAIR
#   else
      thole, pdamp,
#   endif
      0);
}
#elif TINKER9_POLPAIR == 2
// ck.py Version 2.0.2
template <class ETYP>
__global__
void ufield_cu1(int n, TINKER_IMAGE_PARAMS, real off,
   const unsigned* restrict uinfo, int nexclude,
   const int (*restrict exclude)[2], const real* restrict exclude_scale,
   const real* restrict x, const real* restrict y, const real* restrict z,
   const Spatial::SortedAtom* restrict sorted, int nakpl,
   const int* restrict iakpl, int niak, const int* restrict iak,
   const int* restrict lst, const real (*restrict uind)[3],
   const real (*restrict uinp)[3], real (*restrict field)[3],
   real (*restrict fieldp)[3], real aewald)
{
   using d::jpolar;
   using d::njpolar;
   using d::pdamp;
   using d::thlval;

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
   __shared__ int jpi[BLOCK_DIM];
   real ukdx;
   real ukdy;
   real ukdz;
   real ukpx;
   real ukpy;
   real ukpz;
   real pdk;
   int jpk;

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
      jpi[klane] = jpolar[i];
      ukdx = uind[k][0];
      ukdy = uind[k][1];
      ukdz = uind[k][2];
      ukpx = uinp[k][0];
      ukpy = uinp[k][1];
      ukpz = uinp[k][2];
      pdk = pdamp[k];
      jpk = jpolar[k];

      constexpr bool incl = true;
      real xr = xk - xi;
      real yr = yk - yi;
      real zr = zk - zi;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         real pga = thlval[njpolar * jpi[klane] + jpk];
         pair_ufield_v2<ETYP>(r2, xr, yr, zr, scalea, aewald, uidx[klane],
            uidy[klane], uidz[klane], uipx[klane], uipy[klane], uipz[klane],
            pdi[klane], pga, ukdx, ukdy, ukdz, ukpx, ukpy, ukpz, pdk, pga, fidx,
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

      uidx[threadIdx.x] = uind[i][0];
      uidy[threadIdx.x] = uind[i][1];
      uidz[threadIdx.x] = uind[i][2];
      uipx[threadIdx.x] = uinp[i][0];
      uipy[threadIdx.x] = uinp[i][1];
      uipz[threadIdx.x] = uinp[i][2];
      pdi[threadIdx.x] = pdamp[i];
      jpi[threadIdx.x] = jpolar[i];
      ukdx = uind[k][0];
      ukdy = uind[k][1];
      ukdz = uind[k][2];
      ukpx = uinp[k][0];
      ukpy = uinp[k][1];
      ukpz = uinp[k][2];
      pdk = pdamp[k];
      jpk = jpolar[k];
      __syncwarp();

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
            real pga = thlval[njpolar * jpi[klane] + jpk];
            pair_ufield_v2<ETYP>(r2, xr, yr, zr, scalea, aewald, uidx[klane],
               uidy[klane], uidz[klane], uipx[klane], uipy[klane], uipz[klane],
               pdi[klane], pga, ukdx, ukdy, ukdz, ukpx, ukpy, ukpz, pdk, pga,
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
   __syncwarp();

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
      jpi[threadIdx.x] = jpolar[i];
      ukdx = uind[k][0];
      ukdy = uind[k][1];
      ukdz = uind[k][2];
      ukpx = uinp[k][0];
      ukpy = uinp[k][1];
      ukpz = uinp[k][2];
      pdk = pdamp[k];
      jpk = jpolar[k];
      __syncwarp();

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
            real pga = thlval[njpolar * jpi[klane] + jpk];
            pair_ufield_v2<ETYP>(r2, xr, yr, zr, scalea, aewald, uidx[klane],
               uidy[klane], uidz[klane], uipx[klane], uipy[klane], uipz[klane],
               pdi[klane], pga, ukdx, ukdy, ukdz, ukpx, ukpy, ukpz, pdk, pga,
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
   __syncwarp();
}

void ufieldEwaldReal_cu(const real (*uind)[3], const real (*uinp)[3],
   real (*field)[3], real (*fieldp)[3])
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
      st.n, TINKER_IMAGE_ARGS, off, st.si4.bit0, nuexclude, uexclude,
      uexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl, st.niak,
      st.iak, st.lst, uind, uinp, field, fieldp, aewald);
}

void ufieldNonEwald_cu(const real (*uind)[3], const real (*uinp)[3],
   real (*field)[3], real (*fieldp)[3])
{
   const auto& st = *mspatial_v2_unit;
   const real off = switchOff(Switch::MPOLE);

   darray::zero(g::q0, n, field, fieldp);
   int ngrid = gpuGridSize(BLOCK_DIM);
   ngrid *= BLOCK_DIM;
   int nparallel = std::max(st.niak, st.nakpl) * WARP_SIZE;
   nparallel = std::max(nparallel, ngrid);
   launch_k1s(g::s0, nparallel, ufield_cu1<NON_EWALD>, //
      st.n, TINKER_IMAGE_ARGS, off, st.si4.bit0, nuexclude, uexclude,
      uexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl, st.niak,
      st.iak, st.lst, uind, uinp, field, fieldp, 0);
}
#endif
}
