#include "add.h"
#include "empole_chgpen.h"
#include "epolar.h"
#include "epolar_chgpen.h"
#include "glob.spatial.h"
#include "image.h"
#include "induce_donly.h"
#include "launch.h"
#include "mdpq.h"
#include "seq_damp_chgpen.h"
#include "seq_triangle.h"
#include "switch.h"
#include "tool/gpu_card.h"


namespace tinker {
__global__
void sparse_precond_cu3(const real (*restrict rsd)[3], real (*restrict zrsd)[3],
                        const real* restrict polarity, int n, real udiag)
{
   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
        i += blockDim.x * gridDim.x) {
      real poli = udiag * polarity[i];
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
         zrsd[i][j] = poli * rsd[i][j];
         // real test = zrsd[i][j];
         // printf("udiag zrsd %14.6e\n", test);
      }
   }
}


// __launch_bounds__(BLOCK_DIM) __global__
// void sparse_precond_cu4(const real (*restrict rsd)[3], real (*restrict
// zrsd)[3],
//                         real *restrict palpha,
//                         const real* restrict polarity, TINKER_IMAGE_PARAMS,
//                         real cutbuf2, int n,
//                         const Spatial::SortedAtom* restrict sorted, int niak,
//                         const int* restrict iak, const int* restrict lst)
// {
//    const int iwarp = (threadIdx.x + blockIdx.x * blockDim.x) / WARP_SIZE;
//    const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
//    const int ilane = threadIdx.x & (WARP_SIZE - 1);


//    struct Data
//    {
//       real3 fkd;
//       real3 rk, ukd;
//       real  alpha, polk;
//    };
//    __shared__ Data data[BLOCK_DIM];


//    for (int iw = iwarp; iw < niak; iw += nwarp) {
//       real3 fid = make_real3(0, 0, 0);
//       int atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
//       real3 ri = make_real3(sorted[atomi].x, sorted[atomi].y,
//       sorted[atomi].z); int i = sorted[atomi].unsorted; real3 uid =
//       make_real3(rsd[i][0], rsd[i][1], rsd[i][2]); real3 zd =
//       make_real3(zrsd[i][0], zrsd[i][1], zrsd[i][2]); real alphai =
//       palpha[i]; real poli = polarity[i];


//       data[threadIdx.x].fkd = make_real3(0, 0, 0);
//       int shatomk = lst[iw * WARP_SIZE + ilane];
//       data[threadIdx.x].rk =
//          make_real3(sorted[shatomk].x, sorted[shatomk].y, sorted[shatomk].z);
//       int shk = sorted[shatomk].unsorted;
//       data[threadIdx.x].ukd = make_real3(rsd[shk][0], rsd[shk][1],
//       rsd[shk][2]); data[threadIdx.x].alpha = palpha[shk];
//       data[threadIdx.x].polk = polarity[shk];


//       for (int j = 0; j < WARP_SIZE; ++j) {
//          int srclane = (ilane + j) & (WARP_SIZE - 1);
//          int klane = srclane + threadIdx.x - ilane;
//          int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
//          real alphak = data[klane].alpha;
//          real3 dr = data[klane].rk - ri;


//          real r2 = image2(dr.x, dr.y, dr.z);
//          if (atomi < atomk && r2 <= cutbuf2) {
//             real r = REAL_SQRT(r2);
//             real dmpik[3];
//             damp_mut(dmpik,r,alphai,alphak);

//             real polik = poli * data[klane].polk;
//             real rr3 = dmpik[1] * polik * REAL_RECIP(r * r2);
//             real rr5 = 3 * dmpik[2] * polik * REAL_RECIP(r * r2 * r2);


//             real c;
//             c = rr5 * dot3(dr, data[klane].ukd);
//             fid += c * dr - rr3 * data[klane].ukd;

//             c = rr5 * dot3(dr, uid);
//             data[klane].fkd += c * dr - rr3 * uid;

//             //printf("1 %16.8e %16.8e %16.8e\n", fid.x, fid.y, fid.z);
//             // printf("2 %16.8e %16.8e %16.8e\n", zd.x, zd.y, zd.z);

//          } // end if (include)
//       }

//       // printf("%16.8e %16.8e %16.8e\n", aa, ab, ac);

//       atomic_add(fid.x, &zrsd[i][0]);
//       atomic_add(fid.y, &zrsd[i][1]);
//       atomic_add(fid.z, &zrsd[i][2]);
//       atomic_add(data[threadIdx.x].fkd.x, &zrsd[shk][0]);
//       atomic_add(data[threadIdx.x].fkd.y, &zrsd[shk][1]);
//       atomic_add(data[threadIdx.x].fkd.z, &zrsd[shk][2]);
//    } // end for (iw)
// }


// __global__
// void sparse_precond_cu5(const real (*restrict rsd)[3], real (*restrict
// zrsd)[3],
//                         real *restrict palpha,
//                         const real* restrict polarity, TINKER_IMAGE_PARAMS,
//                         real cutbuf2, const real* restrict x,
//                         const real* restrict y, const real* restrict z,
//                         int nwexclude, const int (*restrict wexclude)[2],
//                         const real* restrict wexclude_scale)
// {
//    for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < nwexclude;
//         ii += blockDim.x * gridDim.x) {
//       int i = wexclude[ii][0];
//       int k = wexclude[ii][1];
//       real wscale = wexclude_scale[ii];


//       real xi = x[i];
//       real yi = y[i];
//       real zi = z[i];
//       real alphai = palpha[i];
//       real poli = polarity[i];

//       real alphak = palpha[k];
//       real xr = x[k] - xi;
//       real yr = y[k] - yi;
//       real zr = z[k] - zi;
//       real r2 = image2(xr, yr, zr);
//       if (r2 <= cutbuf2) {
//          real r = REAL_SQRT(r2);
//          real dmpik[3];
//          damp_mut(dmpik,r,alphai,alphak);
//          real scale3 = wscale * dmpik[1];
//          real scale5 = wscale * dmpik[2];

//          real polik = poli * polarity[k];
//          real rr3 = scale3 * polik * REAL_RECIP(r * r2);
//          real rr5 = 3 * scale5 * polik * REAL_RECIP(r * r2 * r2);


//          real c;
//          real3 dr = make_real3(xr, yr, zr);
//          real3 uid = make_real3(rsd[i][0], rsd[i][1], rsd[i][2]);
//          real3 ukd = make_real3(rsd[k][0], rsd[k][1], rsd[k][2]);


//          c = rr5 * dot3(dr, ukd);
//          real3 fid = c * dr - rr3 * ukd;
//          c = rr5 * dot3(dr, uid);
//          real3 fkd = c * dr - rr3 * uid;

//          atomic_add(fid.x, &zrsd[i][0]);
//          atomic_add(fid.y, &zrsd[i][1]);
//          atomic_add(fid.z, &zrsd[i][2]);
//          atomic_add(fkd.x, &zrsd[k][0]);
//          atomic_add(fkd.y, &zrsd[k][1]);
//          atomic_add(fkd.z, &zrsd[k][2]);

//       }
//    }
// }


__launch_bounds__(BLOCK_DIM) __global__
void sparse_precond_cu4(int n, TINKER_IMAGE_PARAMS, real off,
                        const unsigned* restrict uinfo, int nexclude,
                        const int (*restrict exclude)[2],
                        const real* restrict exclude_scale,
                        const real* restrict x, const real* restrict y,
                        const real* restrict z,
                        const Spatial::SortedAtom* restrict sorted, int nakpl,
                        const int* restrict iakpl, int niak,
                        const int* restrict iak, const int* restrict lst,
                        const real (*restrict rsd)[3], real (*restrict zrsd)[3],
                        const real* restrict palpha,
                        const real* restrict polarity)
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
   real fkdx;
   real fkdy;
   real fkdz;
   __shared__ real shuidx[BLOCK_DIM];
   __shared__ real shuidy[BLOCK_DIM];
   __shared__ real shuidz[BLOCK_DIM];
   __shared__ real shalphai[BLOCK_DIM];
   __shared__ real shpoli[BLOCK_DIM];
   real ukdx;
   real ukdy;
   real ukdz;
   real alphak;
   real polk;


   //* /
   // exclude
   for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
      shfidx[threadIdx.x] = 0;
      shfidy[threadIdx.x] = 0;
      shfidz[threadIdx.x] = 0;
      fkdx = 0;
      fkdy = 0;
      fkdz = 0;


      int shi = exclude[ii][0];
      int k = exclude[ii][1];
      real scalea = exclude_scale[ii];


      real xi = x[shi];
      real yi = y[shi];
      real zi = z[shi];
      xk = x[k];
      yk = y[k];
      zk = z[k];
      real uidx = rsd[shi][0];
      real uidy = rsd[shi][1];
      real uidz = rsd[shi][2];
      real alphai = palpha[shi];
      real poli = polarity[shi];
      ukdx = rsd[k][0];
      ukdy = rsd[k][1];
      ukdz = rsd[k][2];
      alphak = palpha[k];
      polk = polarity[k];


      constexpr bool incl = true;
      const int srclane = ilane;
      real xr = xk - xi;
      real yr = yk - yi;
      real zr = zk - zi;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         real r = REAL_SQRT(r2);
         real dmpik[3];
         damp_mut(dmpik, r, alphai, alphak);
         real scale3 = scalea * dmpik[1];
         real scale5 = scalea * dmpik[2];

         real polik = poli * polk;
         real rr3 = scale3 * polik * REAL_RECIP(r * r2);
         real rr5 = 3 * scale5 * polik * REAL_RECIP(r * r2 * r2);

         if (scalea != 0)
         printf(" p1 r rr3 rr5 %16.8e%16.8e%16.8e scale %.4lf\n", r,rr3,rr5,scalea);

         real c;
         int klane = srclane + threadIdx.x - ilane;
         c = rr5 * dot3(xr, yr, zr, ukdx, ukdy, ukdz);
         shfidx[klane] += c * xr - rr3 * ukdx;
         shfidy[klane] += c * yr - rr3 * ukdy;
         shfidz[klane] += c * zr - rr3 * ukdz;


         c = rr5 * dot3(xr, yr, zr, uidx, uidy, uidz);
         fkdx += c * xr - rr3 * shuidx[klane];
         fkdy += c * yr - rr3 * shuidy[klane];
         fkdz += c * zr - rr3 * shuidz[klane];
      } // end if (include)


      atomic_add(shfidx[threadIdx.x], &zrsd[shi][0]);
      atomic_add(shfidy[threadIdx.x], &zrsd[shi][1]);
      atomic_add(shfidz[threadIdx.x], &zrsd[shi][2]);
      atomic_add(fkdx, &zrsd[k][0]);
      atomic_add(fkdy, &zrsd[k][1]);
      atomic_add(fkdz, &zrsd[k][2]);
   }
   // */


   //* /
   // block pairs that have scale factors
   for (int iw = iwarp; iw < nakpl; iw += nwarp) {
      shfidx[threadIdx.x] = 0;
      shfidy[threadIdx.x] = 0;
      shfidz[threadIdx.x] = 0;
      fkdx = 0;
      fkdy = 0;
      fkdz = 0;


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


      shuidx[threadIdx.x] = rsd[shi][0];
      shuidy[threadIdx.x] = rsd[shi][1];
      shuidz[threadIdx.x] = rsd[shi][2];
      shalphai[threadIdx.x] = palpha[shi];
      shpoli[threadIdx.x] = polarity[shi];
      ukdx = rsd[k][0];
      ukdy = rsd[k][1];
      ukdz = rsd[k][2];
      alphak = palpha[k];
      polk = polarity[k];


      unsigned int uinfo0 = uinfo[iw * WARP_SIZE + ilane];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int srcmask = 1 << srclane;
         int iid = shiid;
         real xi = shxi[srclane + threadIdx.x - ilane];
         real yi = shyi[srclane + threadIdx.x - ilane];
         real zi = shzi[srclane + threadIdx.x - ilane];
         real uidx = shuidx[srclane + threadIdx.x - ilane];
         real uidy = shuidy[srclane + threadIdx.x - ilane];
         real uidz = shuidz[srclane + threadIdx.x - ilane];
         real alphai = shalphai[srclane + threadIdx.x - ilane];
         real poli = shpoli[srclane + threadIdx.x - ilane];


         bool incl = iid < kid and kid < n;
         incl = incl and (uinfo0 & srcmask) == 0;
         real scalea = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real r = REAL_SQRT(r2);
            real dmpik[3];
            damp_mut(dmpik, r, alphai, alphak);
            real scale3 = scalea * dmpik[1];
            real scale5 = scalea * dmpik[2];

            real polik = poli * polk;
            real rr3 = scale3 * polik * REAL_RECIP(r * r2);
            real rr5 = 3 * scale5 * polik * REAL_RECIP(r * r2 * r2);

            printf(" p2 r rr3 %16.8e%16.8e%16.8e%16.8e%16.8e%16.8e%16.8e \n", r,rr3,scale3,dmpik[1],poli,polk,polik);

            real c;
            int klane = srclane + threadIdx.x - ilane;
            c = rr5 * dot3(xr, yr, zr, ukdx, ukdy, ukdz);
            shfidx[klane] += c * xr - rr3 * ukdx;
            shfidy[klane] += c * yr - rr3 * ukdy;
            shfidz[klane] += c * zr - rr3 * ukdz;


            c = rr5 * dot3(xr, yr, zr, uidx, uidy, uidz);
            fkdx += c * xr - rr3 * shuidx[klane];
            fkdy += c * yr - rr3 * shuidy[klane];
            fkdz += c * zr - rr3 * shuidz[klane];
         } // end if (include)


         shiid = __shfl_sync(ALL_LANES, shiid, ilane + 1);
      }


      atomic_add(shfidx[threadIdx.x], &zrsd[shi][0]);
      atomic_add(shfidy[threadIdx.x], &zrsd[shi][1]);
      atomic_add(shfidz[threadIdx.x], &zrsd[shi][2]);
      atomic_add(fkdx, &zrsd[k][0]);
      atomic_add(fkdy, &zrsd[k][1]);
      atomic_add(fkdz, &zrsd[k][2]);
   }
   // */


   //* /
   // block-atoms
   for (int iw = iwarp; iw < niak; iw += nwarp) {
      shfidx[threadIdx.x] = 0;
      shfidy[threadIdx.x] = 0;
      shfidz[threadIdx.x] = 0;
      fkdx = 0;
      fkdy = 0;
      fkdz = 0;


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


      shuidx[threadIdx.x] = rsd[shi][0];
      shuidy[threadIdx.x] = rsd[shi][1];
      shuidz[threadIdx.x] = rsd[shi][2];
      shalphai[threadIdx.x] = palpha[shi];
      shpoli[threadIdx.x] = polarity[shi];
      ukdx = rsd[k][0];
      ukdy = rsd[k][1];
      ukdz = rsd[k][2];
      alphak = palpha[k];
      polk = polarity[k];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         real xi = shxi[srclane + threadIdx.x - ilane];
         real yi = shyi[srclane + threadIdx.x - ilane];
         real zi = shzi[srclane + threadIdx.x - ilane];
         real uidx = shuidx[srclane + threadIdx.x - ilane];
         real uidy = shuidy[srclane + threadIdx.x - ilane];
         real uidz = shuidz[srclane + threadIdx.x - ilane];
         real alphai = shalphai[srclane + threadIdx.x - ilane];
         real poli = shpoli[srclane + threadIdx.x - ilane];


         bool incl = atomk > 0;
         real scalea = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            real r = REAL_SQRT(r2);
            real dmpik[3];
            damp_mut(dmpik, r, alphai, alphak);
            real scale3 = scalea * dmpik[1];
            real scale5 = scalea * dmpik[2];

            real polik = poli * polk;
            real rr3 = scale3 * polik * REAL_RECIP(r * r2);
            real rr5 = 3 * scale5 * polik * REAL_RECIP(r * r2 * r2);

            printf(" p3 r rr3 rr5 %16.8e%16.8e%16.8e scale %.4lf\n", r,rr3,rr5,scalea);

            real c;
            int klane = srclane + threadIdx.x - ilane;
            c = rr5 * dot3(xr, yr, zr, ukdx, ukdy, ukdz);
            shfidx[klane] += c * xr - rr3 * ukdx;
            shfidy[klane] += c * yr - rr3 * ukdy;
            shfidz[klane] += c * zr - rr3 * ukdz;


            c = rr5 * dot3(xr, yr, zr, uidx, uidy, uidz);
            fkdx += c * xr - rr3 * shuidx[klane];
            fkdy += c * yr - rr3 * shuidy[klane];
            fkdz += c * zr - rr3 * shuidz[klane];
         } // end if (include)
      }


      atomic_add(shfidx[threadIdx.x], &zrsd[shi][0]);
      atomic_add(shfidy[threadIdx.x], &zrsd[shi][1]);
      atomic_add(shfidz[threadIdx.x], &zrsd[shi][2]);
      atomic_add(fkdx, &zrsd[k][0]);
      atomic_add(fkdy, &zrsd[k][1]);
      atomic_add(fkdz, &zrsd[k][2]);
   }
   // */
} // generated by ComplexKernelBuilder (ck.py) 1.3


void sparse_precond_apply_cu2(const real (*rsd)[3], real (*zrsd)[3])
{
   // const auto& st = *uspatial_unit;
   // const real off = switch_off(switch_usolve);
   // const real cutbuf2 = (off + st.buffer) * (off + st.buffer);


   const auto& st = *uspatial_v2_unit;
   const real off = switch_off(switch_usolve) + st.buffer;


   launch_k1s(nonblk, n, sparse_precond_cu3, //
              rsd, zrsd, polarity, n, udiag);
   // if (st.niak > 0)
   //    launch_k1s(nonblk, WARP_SIZE * st.niak, sparse_precond_cu4, //
   //               rsd, zrsd, palpha, polarity, TINKER_IMAGE_ARGS,
   //               cutbuf2, //
   //               n, st.sorted, st.niak, st.iak, st.lst);

   // if (nwexclude > 0)
   //    launch_k1s(nonblk, nwexclude, sparse_precond_cu5, //
   //               rsd, zrsd, palpha, polarity, TINKER_IMAGE_ARGS,
   //               cutbuf2, //
   //               x, y, z, nwexclude, wexclude, wexclude_scale);
   int ngrid = get_grid_size(BLOCK_DIM);
   sparse_precond_cu4<<<ngrid, BLOCK_DIM, 0, nonblk>>>(
      st.n, TINKER_IMAGE_ARGS, off, st.si1.bit0, nuexclude, uexclude,
      uexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl, st.niak,
      st.iak, st.lst, rsd, zrsd, palpha, polarity);
}
}
