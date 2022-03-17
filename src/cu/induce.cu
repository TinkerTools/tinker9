#include "add.h"
#include "epolar.h"
#include "glob.spatial.h"
#include "image.h"
#include "induce.h"
#include "launch.h"
#include "md.h"
#include "seq_damp.h"
#include "seq_triangle.h"
#include "spatial2.h"
#include "switch.h"
#include "tool/gpu_card.h"

namespace tinker {
__global__
void sparse_precond_cu0(const real (*restrict rsd)[3], const real (*restrict rsdp)[3],
   real (*restrict zrsd)[3], real (*restrict zrsdp)[3], const real* restrict polarity, int n,
   real udiag)
{
   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n; i += blockDim.x * gridDim.x) {
      real poli = udiag * polarity[i];
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
         zrsd[i][j] = poli * rsd[i][j];
         zrsdp[i][j] = poli * rsdp[i][j];
      }
   }
}

// ck.py Version 2.0.2
__global__
void sparse_precond_cu1(int n, TINKER_IMAGE_PARAMS, real off, const unsigned* restrict uinfo,
   int nexclude, const int (*restrict exclude)[2], const real* restrict exclude_scale,
   const real* restrict x, const real* restrict y, const real* restrict z,
   const Spatial::SortedAtom* restrict sorted, int nakpl, const int* restrict iakpl, int niak,
   const int* restrict iak, const int* restrict lst, const real (*restrict rsd)[3],
   const real (*restrict rsdp)[3], real (*restrict zrsd)[3], real (*restrict zrsdp)[3],
   const real* restrict pdamp, const real* restrict thole, const real* restrict polarity)
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
   __shared__ real poli[BLOCK_DIM];
   real ukdx;
   real ukdy;
   real ukdz;
   real ukpx;
   real ukpy;
   real ukpz;
   real pdk;
   real ptk;
   real polk;

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
      uidx[klane] = rsd[i][0];
      uidy[klane] = rsd[i][1];
      uidz[klane] = rsd[i][2];
      uipx[klane] = rsdp[i][0];
      uipy[klane] = rsdp[i][1];
      uipz[klane] = rsdp[i][2];
      pdi[klane] = pdamp[i];
      pti[klane] = thole[i];
      poli[klane] = polarity[i];
      ukdx = rsd[k][0];
      ukdy = rsd[k][1];
      ukdz = rsd[k][2];
      ukpx = rsdp[k][0];
      ukpy = rsdp[k][1];
      ukpz = rsdp[k][2];
      pdk = pdamp[k];
      ptk = thole[k];
      polk = polarity[k];

      constexpr bool incl = true;
      real xr = xk - xi;
      real yr = yk - yi;
      real zr = zk - zi;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         real r = REAL_SQRT(r2);
         real scale3, scale5;
         damp_thole2(r, pdi[klane], pti[klane], pdk, ptk, scale3, scale5);
         scale3 *= scalea;
         scale5 *= scalea;
         real polik = poli[klane] * polk;
         real rr3 = scale3 * polik * REAL_RECIP(r * r2);
         real rr5 = 3 * scale5 * polik * REAL_RECIP(r * r2 * r2);

         real c;
         c = rr5 * dot3(xr, yr, zr, ukdx, ukdy, ukdz);
         fidx += c * xr - rr3 * ukdx;
         fidz += c * zr - rr3 * ukdz;
         fidy += c * yr - rr3 * ukdy;

         c = rr5 * dot3(xr, yr, zr, ukpx, ukpy, ukpz);
         fipx += c * xr - rr3 * ukpx;
         fipy += c * yr - rr3 * ukpy;
         fipz += c * zr - rr3 * ukpz;

         c = rr5 * dot3(xr, yr, zr, uidx[klane], uidy[klane], uidz[klane]);
         fkdx += c * xr - rr3 * uidx[klane];
         fkdy += c * yr - rr3 * uidy[klane];
         fkdz += c * zr - rr3 * uidz[klane];

         c = rr5 * dot3(xr, yr, zr, uipx[klane], uipy[klane], uipz[klane]);
         fkpx += c * xr - rr3 * uipx[klane];
         fkpy += c * yr - rr3 * uipy[klane];
         fkpz += c * zr - rr3 * uipz[klane];
      } // end if (include)

      atomic_add(fidx, &zrsd[i][0]);
      atomic_add(fidy, &zrsd[i][1]);
      atomic_add(fidz, &zrsd[i][2]);
      atomic_add(fipx, &zrsdp[i][0]);
      atomic_add(fipy, &zrsdp[i][1]);
      atomic_add(fipz, &zrsdp[i][2]);
      atomic_add(fkdx, &zrsd[k][0]);
      atomic_add(fkdy, &zrsd[k][1]);
      atomic_add(fkdz, &zrsd[k][2]);
      atomic_add(fkpx, &zrsdp[k][0]);
      atomic_add(fkpy, &zrsdp[k][1]);
      atomic_add(fkpz, &zrsdp[k][2]);
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

      uidx[threadIdx.x] = rsd[i][0];
      uidy[threadIdx.x] = rsd[i][1];
      uidz[threadIdx.x] = rsd[i][2];
      uipx[threadIdx.x] = rsdp[i][0];
      uipy[threadIdx.x] = rsdp[i][1];
      uipz[threadIdx.x] = rsdp[i][2];
      pdi[threadIdx.x] = pdamp[i];
      pti[threadIdx.x] = thole[i];
      poli[threadIdx.x] = polarity[i];
      ukdx = rsd[k][0];
      ukdy = rsd[k][1];
      ukdz = rsd[k][2];
      ukpx = rsdp[k][0];
      ukpy = rsdp[k][1];
      ukpz = rsdp[k][2];
      pdk = pdamp[k];
      ptk = thole[k];
      polk = polarity[k];

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
            real r = REAL_SQRT(r2);
            real scale3, scale5;
            damp_thole2(r, pdi[klane], pti[klane], pdk, ptk, scale3, scale5);
            scale3 *= scalea;
            scale5 *= scalea;
            real polik = poli[klane] * polk;
            real rr3 = scale3 * polik * REAL_RECIP(r * r2);
            real rr5 = 3 * scale5 * polik * REAL_RECIP(r * r2 * r2);

            real c;
            c = rr5 * dot3(xr, yr, zr, ukdx, ukdy, ukdz);
            fidx += c * xr - rr3 * ukdx;
            fidz += c * zr - rr3 * ukdz;
            fidy += c * yr - rr3 * ukdy;

            c = rr5 * dot3(xr, yr, zr, ukpx, ukpy, ukpz);
            fipx += c * xr - rr3 * ukpx;
            fipy += c * yr - rr3 * ukpy;
            fipz += c * zr - rr3 * ukpz;

            c = rr5 * dot3(xr, yr, zr, uidx[klane], uidy[klane], uidz[klane]);
            fkdx += c * xr - rr3 * uidx[klane];
            fkdy += c * yr - rr3 * uidy[klane];
            fkdz += c * zr - rr3 * uidz[klane];

            c = rr5 * dot3(xr, yr, zr, uipx[klane], uipy[klane], uipz[klane]);
            fkpx += c * xr - rr3 * uipx[klane];
            fkpy += c * yr - rr3 * uipy[klane];
            fkpz += c * zr - rr3 * uipz[klane];
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

      atomic_add(fidx, &zrsd[i][0]);
      atomic_add(fidy, &zrsd[i][1]);
      atomic_add(fidz, &zrsd[i][2]);
      atomic_add(fipx, &zrsdp[i][0]);
      atomic_add(fipy, &zrsdp[i][1]);
      atomic_add(fipz, &zrsdp[i][2]);
      atomic_add(fkdx, &zrsd[k][0]);
      atomic_add(fkdy, &zrsd[k][1]);
      atomic_add(fkdz, &zrsd[k][2]);
      atomic_add(fkpx, &zrsdp[k][0]);
      atomic_add(fkpy, &zrsdp[k][1]);
      atomic_add(fkpz, &zrsdp[k][2]);
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

      uidx[threadIdx.x] = rsd[i][0];
      uidy[threadIdx.x] = rsd[i][1];
      uidz[threadIdx.x] = rsd[i][2];
      uipx[threadIdx.x] = rsdp[i][0];
      uipy[threadIdx.x] = rsdp[i][1];
      uipz[threadIdx.x] = rsdp[i][2];
      pdi[threadIdx.x] = pdamp[i];
      pti[threadIdx.x] = thole[i];
      poli[threadIdx.x] = polarity[i];
      ukdx = rsd[k][0];
      ukdy = rsd[k][1];
      ukdz = rsd[k][2];
      ukpx = rsdp[k][0];
      ukpy = rsdp[k][1];
      ukpz = rsdp[k][2];
      pdk = pdamp[k];
      ptk = thole[k];
      polk = polarity[k];

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
            real r = REAL_SQRT(r2);
            real scale3, scale5;
            damp_thole2(r, pdi[klane], pti[klane], pdk, ptk, scale3, scale5);
            scale3 *= scalea;
            scale5 *= scalea;
            real polik = poli[klane] * polk;
            real rr3 = scale3 * polik * REAL_RECIP(r * r2);
            real rr5 = 3 * scale5 * polik * REAL_RECIP(r * r2 * r2);

            real c;
            c = rr5 * dot3(xr, yr, zr, ukdx, ukdy, ukdz);
            fidx += c * xr - rr3 * ukdx;
            fidz += c * zr - rr3 * ukdz;
            fidy += c * yr - rr3 * ukdy;

            c = rr5 * dot3(xr, yr, zr, ukpx, ukpy, ukpz);
            fipx += c * xr - rr3 * ukpx;
            fipy += c * yr - rr3 * ukpy;
            fipz += c * zr - rr3 * ukpz;

            c = rr5 * dot3(xr, yr, zr, uidx[klane], uidy[klane], uidz[klane]);
            fkdx += c * xr - rr3 * uidx[klane];
            fkdy += c * yr - rr3 * uidy[klane];
            fkdz += c * zr - rr3 * uidz[klane];

            c = rr5 * dot3(xr, yr, zr, uipx[klane], uipy[klane], uipz[klane]);
            fkpx += c * xr - rr3 * uipx[klane];
            fkpy += c * yr - rr3 * uipy[klane];
            fkpz += c * zr - rr3 * uipz[klane];
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

      atomic_add(fidx, &zrsd[i][0]);
      atomic_add(fidy, &zrsd[i][1]);
      atomic_add(fidz, &zrsd[i][2]);
      atomic_add(fipx, &zrsdp[i][0]);
      atomic_add(fipy, &zrsdp[i][1]);
      atomic_add(fipz, &zrsdp[i][2]);
      atomic_add(fkdx, &zrsd[k][0]);
      atomic_add(fkdy, &zrsd[k][1]);
      atomic_add(fkdz, &zrsd[k][2]);
      atomic_add(fkpx, &zrsdp[k][0]);
      atomic_add(fkpy, &zrsdp[k][1]);
      atomic_add(fkpz, &zrsdp[k][2]);
   }
}

void sparse_precond_apply_cu(
   const real (*rsd)[3], const real (*rsdp)[3], real (*zrsd)[3], real (*zrsdp)[3])
{
   const auto& st = *uspatial_v2_unit;
   real off = switch_off(switch_usolve);
   off = off + st.buffer;

   launch_k1s(g::s0, n, sparse_precond_cu0, //
      rsd, rsdp, zrsd, zrsdp, polarity, n, udiag);
   int ngrid = gpuGridSize(BLOCK_DIM);
   ngrid *= BLOCK_DIM;
   int nparallel = std::max(st.niak, st.nakpl) * WARP_SIZE;
   nparallel = std::max(nparallel, ngrid);
   launch_k1s(g::s0, nparallel, sparse_precond_cu1, //
      st.n, TINKER_IMAGE_ARGS, off, st.si1.bit0, nuexclude, uexclude, uexclude_scale, st.x, st.y,
      st.z, st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst, rsd, rsdp, zrsd, zrsdp, pdamp,
      thole, polarity);
}
}
