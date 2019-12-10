#include "add.cuh"
#include "induce.h"
#include "launch.cuh"
#include "spatial.h"


TINKER_NAMESPACE_BEGIN
__global__
void sparse_precond_cu0(const real (*restrict rsd)[3],
                        const real (*restrict rsdp)[3],
                        real (*restrict zrsd)[3], real (*restrict zrsdp)[3],
                        const real* restrict polarity, int n, real udiag)
{
   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
        i += blockDim.x * gridDim.x) {
      real poli = udiag * polarity[i];
      for (int j = 0; j < 3; ++j) {
         zrsd[i][j] = poli * rsd[i][j];
         zrsdp[i][j] = poli * rsdp[i][j];
      }
   }
}


__global__
void sparse_precond_cu1(const real (*restrict rsd)[3],
                        const real (*restrict rsdp)[3],
                        real (*restrict zrsd)[3], real (*restrict zrsdp)[3],
                        const real* restrict pdamp, const real* restrict thole,
                        const real* restrict polarity, const Box* restrict box,
                        real cutbuf2, int n,
                        const Spatial::SortedAtom* restrict sorted, int niak,
                        const int* restrict iak, const int* restrict lst)
{
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);


   real gxi, gyi, gzi, txi, tyi, tzi;
   __shared__ real gxk[BLOCK_DIM], gyk[BLOCK_DIM], gzk[BLOCK_DIM],
      txk[BLOCK_DIM], tyk[BLOCK_DIM], tzk[BLOCK_DIM];


   for (int iw = iwarp; iw < niak; iw += nwarp) {
      gxi = 0;
      gyi = 0;
      gzi = 0;
      txi = 0;
      tyi = 0;
      tzi = 0;
      gxk[threadIdx.x] = 0;
      gyk[threadIdx.x] = 0;
      gzk[threadIdx.x] = 0;
      txk[threadIdx.x] = 0;
      tyk[threadIdx.x] = 0;
      tzk[threadIdx.x] = 0;


      int atomi;
      atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
      real xi = sorted[atomi].x;
      real yi = sorted[atomi].y;
      real zi = sorted[atomi].z;
      int i = sorted[atomi].unsorted;
      real pdi = pdamp[i];
      real pti = thole[i];
      real poli = polarity[i];


      int shatomk;
      shatomk = lst[iw * WARP_SIZE + ilane];
      real shx = sorted[shatomk].x;
      real shy = sorted[shatomk].y;
      real shz = sorted[shatomk].z;
      int shk = sorted[shatomk].unsorted;
      real shpdk = pdamp[shk];
      real shptk = thole[shk];
      real shpolk = polarity[shk];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
         real xr = __shfl_sync(ALL_LANES, shx, srclane) - xi;
         real yr = __shfl_sync(ALL_LANES, shy, srclane) - yi;
         real zr = __shfl_sync(ALL_LANES, shz, srclane) - zi;
         int k = __shfl_sync(ALL_LANES, shk, srclane);
         real pdk = __shfl_sync(ALL_LANES, shpdk, srclane);
         real ptk = __shfl_sync(ALL_LANES, shptk, srclane);
         real polk = __shfl_sync(ALL_LANES, shpolk, srclane);


         real m0 = 0;
         real m1 = 0;
         real m2 = 0;
         real m3 = 0;
         real m4 = 0;
         real m5 = 0;


         image(xr, yr, zr, box);
         real r2 = xr * xr + yr * yr + zr * zr;
         if (atomi < atomk && r2 <= cutbuf2) {
            real r = REAL_SQRT(r2);
            real scale3, scale5;
            damp_thole2(r, pdi, pti, pdk, ptk, scale3, scale5);
            real polik = poli * polk;
            real rr3 = scale3 * polik * REAL_RECIP(r * r2);
            real rr5 = 3 * scale5 * polik * REAL_RECIP(r * r2 * r2);


            m0 = rr5 * xr * xr - rr3;
            m1 = rr5 * xr * yr;
            m2 = rr5 * xr * zr;
            m3 = rr5 * yr * yr - rr3;
            m4 = rr5 * yr * zr;
            m5 = rr5 * zr * zr - rr3;
         } // end if (include)


         gxi += m0 * rsd[k][0] + m1 * rsd[k][1] + m2 * rsd[k][2];
         gyi += m1 * rsd[k][0] + m3 * rsd[k][1] + m4 * rsd[k][2];
         gzi += m2 * rsd[k][0] + m4 * rsd[k][1] + m5 * rsd[k][2];
         gxk[srclane + (threadIdx.x - ilane)] +=
            m0 * rsd[i][0] + m1 * rsd[i][1] + m2 * rsd[i][2];
         gyk[srclane + (threadIdx.x - ilane)] +=
            m1 * rsd[i][0] + m3 * rsd[i][1] + m4 * rsd[i][2];
         gzk[srclane + (threadIdx.x - ilane)] +=
            m2 * rsd[i][0] + m4 * rsd[i][1] + m5 * rsd[i][2];
         txi += m0 * rsdp[k][0] + m1 * rsdp[k][1] + m2 * rsdp[k][2];
         tyi += m1 * rsdp[k][0] + m3 * rsdp[k][1] + m4 * rsdp[k][2];
         tzi += m2 * rsdp[k][0] + m4 * rsdp[k][1] + m5 * rsdp[k][2];
         txk[srclane + (threadIdx.x - ilane)] +=
            m0 * rsdp[i][0] + m1 * rsdp[i][1] + m2 * rsdp[i][2];
         tyk[srclane + (threadIdx.x - ilane)] +=
            m1 * rsdp[i][0] + m3 * rsdp[i][1] + m4 * rsdp[i][2];
         tzk[srclane + (threadIdx.x - ilane)] +=
            m2 * rsdp[i][0] + m4 * rsdp[i][1] + m5 * rsdp[i][2];
      } // end for (j)

      atomic_add(gxi, &zrsd[i][0]);
      atomic_add(gyi, &zrsd[i][1]);
      atomic_add(gzi, &zrsd[i][2]);
      atomic_add(txi, &zrsdp[i][0]);
      atomic_add(tyi, &zrsdp[i][1]);
      atomic_add(tzi, &zrsdp[i][2]);
      atomic_add(gxk[threadIdx.x], &zrsd[shk][0]);
      atomic_add(gyk[threadIdx.x], &zrsd[shk][1]);
      atomic_add(gzk[threadIdx.x], &zrsd[shk][2]);
      atomic_add(txk[threadIdx.x], &zrsdp[shk][0]);
      atomic_add(tyk[threadIdx.x], &zrsdp[shk][1]);
      atomic_add(tzk[threadIdx.x], &zrsdp[shk][2]);
   } // end for (iw)
}


__global__
void sparse_precond_cu2(const real (*restrict rsd)[3],
                        const real (*restrict rsdp)[3],
                        real (*restrict zrsd)[3], real (*restrict zrsdp)[3],
                        const real* restrict pdamp, const real* restrict thole,
                        const real* restrict polarity, const Box* restrict box,
                        real cutbuf2, const real* restrict x,
                        const real* restrict y, const real* restrict z,
                        int nuexclude_, const int (*restrict uexclude_)[2],
                        const real* restrict uexclude_scale_)
{
   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < nuexclude_;
        ii += blockDim.x * gridDim.x) {
      int i = uexclude_[ii][0];
      int k = uexclude_[ii][1];
      real uscale = uexclude_scale_[ii];


      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real pdi = pdamp[i];
      real pti = thole[i];
      real poli = polarity[i];


      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;
      image(xr, yr, zr, box);
      real r2 = xr * xr + yr * yr + zr * zr;
      if (r2 <= cutbuf2) {
         real r = REAL_SQRT(r2);
         real scale3, scale5;
         damp_thole2(r, pdi, pti, pdamp[k], thole[k], scale3, scale5);
         scale3 *= uscale;
         scale5 *= uscale;
         real polik = poli * polarity[k];
         real rr3 = scale3 * polik * REAL_RECIP(r * r2);
         real rr5 = 3 * scale5 * polik * REAL_RECIP(r * r2 * r2);


         real m0 = rr5 * xr * xr - rr3;
         real m1 = rr5 * xr * yr;
         real m2 = rr5 * xr * zr;
         real m3 = rr5 * yr * yr - rr3;
         real m4 = rr5 * yr * zr;
         real m5 = rr5 * zr * zr - rr3;
         atomic_add(m0 * rsd[k][0] + m1 * rsd[k][1] + m2 * rsd[k][2],
                    &zrsd[i][0]);
         atomic_add(m1 * rsd[k][0] + m3 * rsd[k][1] + m4 * rsd[k][2],
                    &zrsd[i][1]);
         atomic_add(m2 * rsd[k][0] + m4 * rsd[k][1] + m5 * rsd[k][2],
                    &zrsd[i][2]);
         atomic_add(m0 * rsd[i][0] + m1 * rsd[i][1] + m2 * rsd[i][2],
                    &zrsd[k][0]);
         atomic_add(m1 * rsd[i][0] + m3 * rsd[i][1] + m4 * rsd[i][2],
                    &zrsd[k][1]);
         atomic_add(m2 * rsd[i][0] + m4 * rsd[i][1] + m5 * rsd[i][2],
                    &zrsd[k][2]);
         atomic_add(m0 * rsdp[k][0] + m1 * rsdp[k][1] + m2 * rsdp[k][2],
                    &zrsdp[i][0]);
         atomic_add(m1 * rsdp[k][0] + m3 * rsdp[k][1] + m4 * rsdp[k][2],
                    &zrsdp[i][1]);
         atomic_add(m2 * rsdp[k][0] + m4 * rsdp[k][1] + m5 * rsdp[k][2],
                    &zrsdp[i][2]);
         atomic_add(m0 * rsdp[i][0] + m1 * rsdp[i][1] + m2 * rsdp[i][2],
                    &zrsdp[k][0]);
         atomic_add(m1 * rsdp[i][0] + m3 * rsdp[i][1] + m4 * rsdp[i][2],
                    &zrsdp[k][1]);
         atomic_add(m2 * rsdp[i][0] + m4 * rsdp[i][1] + m5 * rsdp[i][2],
                    &zrsdp[k][2]);
      }
   }
}


void sparse_precond_apply_cu(const real (*rsd)[3], const real (*rsdp)[3],
                             real (*zrsd)[3], real (*zrsdp)[3])
{
   const auto& st = *uspatial_unit;
   const real cutbuf2 = (st.cutoff + st.buffer) * (st.cutoff + st.buffer);

   launch_kernel1(n, sparse_precond_cu0, //
                  rsd, rsdp, zrsd, zrsdp, polarity, n, udiag);
   if (st.niak > 0)
      launch_kernel1(WARP_SIZE * st.niak, sparse_precond_cu1, //
                     rsd, rsdp, zrsd, zrsdp, pdamp, thole, polarity, box,
                     cutbuf2, //
                     n, st.sorted, st.niak, st.iak, st.lst);
   if (nuexclude_ > 0)
      launch_kernel1(nuexclude_, sparse_precond_cu2, //
                     rsd, rsdp, zrsd, zrsdp, pdamp, thole, polarity, box,
                     cutbuf2, //
                     x, y, z, nuexclude_, uexclude_, uexclude_scale_);
}
TINKER_NAMESPACE_END
