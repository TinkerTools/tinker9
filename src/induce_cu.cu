#include "add.h"
#include "epolar.h"
#include "image.h"
#include "induce.h"
#include "launch.h"
#include "mdpq.h"
#include "seq_damp.h"
#include "spatial.h"
#include "switch.h"


namespace tinker {
__global__
void sparse_precond_cu0(const real (*restrict rsd)[3],
                        const real (*restrict rsdp)[3],
                        real (*restrict zrsd)[3], real (*restrict zrsdp)[3],
                        const real* restrict polarity, int n, real udiag)
{
   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
        i += blockDim.x * gridDim.x) {
      real poli = udiag * polarity[i];
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
         zrsd[i][j] = poli * rsd[i][j];
         zrsdp[i][j] = poli * rsdp[i][j];
      }
   }
}


__launch_bounds__(BLOCK_DIM) __global__
void sparse_precond_cu1(const real (*restrict rsd)[3],
                        const real (*restrict rsdp)[3],
                        real (*restrict zrsd)[3], real (*restrict zrsdp)[3],
                        const real* restrict pdamp, const real* restrict thole,
                        const real* restrict polarity, TINKER_IMAGE_PARAMS,
                        real cutbuf2, int n,
                        const Spatial::SortedAtom* restrict sorted, int niak,
                        const int* restrict iak, const int* restrict lst)
{
   const int iwarp = (threadIdx.x + blockIdx.x * blockDim.x) / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);


   struct Data
   {
      real3 fkd, fkp;
      real3 rk, ukd, ukp;
      real pdk, ptk, polk;
   };
   __shared__ Data data[BLOCK_DIM];


   for (int iw = iwarp; iw < niak; iw += nwarp) {
      real3 fid = make_real3(0, 0, 0);
      real3 fip = make_real3(0, 0, 0);
      int atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
      real3 ri = make_real3(sorted[atomi].x, sorted[atomi].y, sorted[atomi].z);
      int i = sorted[atomi].unsorted;
      real3 uid = make_real3(rsd[i][0], rsd[i][1], rsd[i][2]);
      real3 uip = make_real3(rsdp[i][0], rsdp[i][1], rsdp[i][2]);
      real pdi = pdamp[i];
      real pti = thole[i];
      real poli = polarity[i];


      data[threadIdx.x].fkd = make_real3(0, 0, 0);
      data[threadIdx.x].fkp = make_real3(0, 0, 0);
      int shatomk = lst[iw * WARP_SIZE + ilane];
      data[threadIdx.x].rk =
         make_real3(sorted[shatomk].x, sorted[shatomk].y, sorted[shatomk].z);
      int shk = sorted[shatomk].unsorted;
      data[threadIdx.x].ukd = make_real3(rsd[shk][0], rsd[shk][1], rsd[shk][2]);
      data[threadIdx.x].ukp =
         make_real3(rsdp[shk][0], rsdp[shk][1], rsdp[shk][2]);
      data[threadIdx.x].pdk = pdamp[shk];
      data[threadIdx.x].ptk = thole[shk];
      data[threadIdx.x].polk = polarity[shk];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
         real3 dr = data[klane].rk - ri;


         real r2 = image2(dr.x, dr.y, dr.z);
         if (atomi < atomk && r2 <= cutbuf2) {
            real r = REAL_SQRT(r2);
            real scale3, scale5;
            damp_thole2(r, pdi, pti, data[klane].pdk, data[klane].ptk, scale3,
                        scale5);
            real polik = poli * data[klane].polk;
            real rr3 = scale3 * polik * REAL_RECIP(r * r2);
            real rr5 = 3 * scale5 * polik * REAL_RECIP(r * r2 * r2);


            real c;
            c = rr5 * dot3(dr, data[klane].ukd);
            fid += c * dr - rr3 * data[klane].ukd;

            c = rr5 * dot3(dr, data[klane].ukp);
            fip += c * dr - rr3 * data[klane].ukp;

            c = rr5 * dot3(dr, uid);
            data[klane].fkd += c * dr - rr3 * uid;

            c = rr5 * dot3(dr, uip);
            data[klane].fkp += c * dr - rr3 * uip;
         } // end if (include)
      }


      atomic_add(fid.x, &zrsd[i][0]);
      atomic_add(fid.y, &zrsd[i][1]);
      atomic_add(fid.z, &zrsd[i][2]);
      atomic_add(fip.x, &zrsdp[i][0]);
      atomic_add(fip.y, &zrsdp[i][1]);
      atomic_add(fip.z, &zrsdp[i][2]);
      atomic_add(data[threadIdx.x].fkd.x, &zrsd[shk][0]);
      atomic_add(data[threadIdx.x].fkd.y, &zrsd[shk][1]);
      atomic_add(data[threadIdx.x].fkd.z, &zrsd[shk][2]);
      atomic_add(data[threadIdx.x].fkp.x, &zrsdp[shk][0]);
      atomic_add(data[threadIdx.x].fkp.y, &zrsdp[shk][1]);
      atomic_add(data[threadIdx.x].fkp.z, &zrsdp[shk][2]);
   } // end for (iw)
}


__global__
void sparse_precond_cu2(const real (*restrict rsd)[3],
                        const real (*restrict rsdp)[3],
                        real (*restrict zrsd)[3], real (*restrict zrsdp)[3],
                        const real* restrict pdamp, const real* restrict thole,
                        const real* restrict polarity, TINKER_IMAGE_PARAMS,
                        real cutbuf2, const real* restrict x,
                        const real* restrict y, const real* restrict z,
                        int nuexclude, const int (*restrict uexclude)[2],
                        const real* restrict uexclude_scale)
{
   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < nuexclude;
        ii += blockDim.x * gridDim.x) {
      int i = uexclude[ii][0];
      int k = uexclude[ii][1];
      real uscale = uexclude_scale[ii];


      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real pdi = pdamp[i];
      real pti = thole[i];
      real poli = polarity[i];


      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;
      real r2 = image2(xr, yr, zr);
      if (r2 <= cutbuf2) {
         real r = REAL_SQRT(r2);
         real scale3, scale5;
         damp_thole2(r, pdi, pti, pdamp[k], thole[k], scale3, scale5);
         scale3 *= uscale;
         scale5 *= uscale;
         real polik = poli * polarity[k];
         real rr3 = scale3 * polik * REAL_RECIP(r * r2);
         real rr5 = 3 * scale5 * polik * REAL_RECIP(r * r2 * r2);


         real c;
         real3 dr = make_real3(xr, yr, zr);
         real3 uid = make_real3(rsd[i][0], rsd[i][1], rsd[i][2]);
         real3 ukd = make_real3(rsd[k][0], rsd[k][1], rsd[k][2]);
         real3 uip = make_real3(rsdp[i][0], rsdp[i][1], rsdp[i][2]);
         real3 ukp = make_real3(rsdp[k][0], rsdp[k][1], rsdp[k][2]);


         c = rr5 * dot3(dr, ukd);
         real3 fid = c * dr - rr3 * ukd;
         c = rr5 * dot3(dr, ukp);
         real3 fip = c * dr - rr3 * ukp;
         c = rr5 * dot3(dr, uid);
         real3 fkd = c * dr - rr3 * uid;
         c = rr5 * dot3(dr, uip);
         real3 fkp = c * dr - rr3 * uip;


         atomic_add(fid.x, &zrsd[i][0]);
         atomic_add(fid.y, &zrsd[i][1]);
         atomic_add(fid.z, &zrsd[i][2]);
         atomic_add(fip.x, &zrsdp[i][0]);
         atomic_add(fip.y, &zrsdp[i][1]);
         atomic_add(fip.z, &zrsdp[i][2]);
         atomic_add(fkd.x, &zrsd[k][0]);
         atomic_add(fkd.y, &zrsd[k][1]);
         atomic_add(fkd.z, &zrsd[k][2]);
         atomic_add(fkp.x, &zrsdp[k][0]);
         atomic_add(fkp.y, &zrsdp[k][1]);
         atomic_add(fkp.z, &zrsdp[k][2]);
      }
   }
}


void sparse_precond_apply_cu(const real (*rsd)[3], const real (*rsdp)[3],
                             real (*zrsd)[3], real (*zrsdp)[3])
{
   const auto& st = *uspatial_unit;
   const real off = switch_off(switch_usolve);
   const real cutbuf2 = (off + st.buffer) * (off + st.buffer);


   launch_k1s(nonblk, n, sparse_precond_cu0, //
              rsd, rsdp, zrsd, zrsdp, polarity, n, udiag);
   if (st.niak > 0)
      launch_k1s(nonblk, WARP_SIZE * st.niak, sparse_precond_cu1, //
                 rsd, rsdp, zrsd, zrsdp, pdamp, thole, polarity,
                 TINKER_IMAGE_ARGS, cutbuf2, //
                 n, st.sorted, st.niak, st.iak, st.lst);
   if (nuexclude > 0)
      launch_k1s(nonblk, nuexclude, sparse_precond_cu2, //
                 rsd, rsdp, zrsd, zrsdp, pdamp, thole, polarity,
                 TINKER_IMAGE_ARGS, cutbuf2, //
                 x, y, z, nuexclude, uexclude, uexclude_scale);
}
}
