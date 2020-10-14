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
      for (int j = 0; j < 3; ++j)
         zrsd[i][j] = poli * rsd[i][j];
   }
}

__launch_bounds__(BLOCK_DIM) __global__
void sparse_precond_cu4(int n, TINKER_IMAGE_PARAMS, real off,
                        const unsigned* restrict winfo, int nexclude,
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


      int klane = threadIdx.x;
      constexpr bool incl = true;
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

         real c;
         c = rr5 * dot3(xr, yr, zr, ukdx, ukdy, ukdz);
         shfidx[klane] += c * xr - rr3 * ukdx;
         shfidy[klane] += c * yr - rr3 * ukdy;
         shfidz[klane] += c * zr - rr3 * ukdz;


         c = rr5 * dot3(xr, yr, zr, uidx, uidy, uidz);
         fkdx += c * xr - rr3 * uidx;
         fkdy += c * yr - rr3 * uidy;
         fkdz += c * zr - rr3 * uidz;
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


      unsigned int winfo0 = winfo[iw * WARP_SIZE + ilane];


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
         real alphai = shalphai[klane];
         real poli = shpoli[klane];


         bool incl = iid < kid and kid < n;
         incl = incl and (winfo0 & srcmask) == 0;
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

            real c;
            c = rr5 * dot3(xr, yr, zr, ukdx, ukdy, ukdz);
            shfidx[klane] += c * xr - rr3 * ukdx;
            shfidy[klane] += c * yr - rr3 * ukdy;
            shfidz[klane] += c * zr - rr3 * ukdz;


            c = rr5 * dot3(xr, yr, zr, uidx, uidy, uidz);
            fkdx += c * xr - rr3 * uidx;
            fkdy += c * yr - rr3 * uidy;
            fkdz += c * zr - rr3 * uidz;
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
         int klane = srclane + threadIdx.x - ilane;
         real xi = shxi[klane];
         real yi = shyi[klane];
         real zi = shzi[klane];
         real uidx = shuidx[klane];
         real uidy = shuidy[klane];
         real uidz = shuidz[klane];
         real alphai = shalphai[klane];
         real poli = shpoli[klane];


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

            real c;
            c = rr5 * dot3(xr, yr, zr, ukdx, ukdy, ukdz);
            shfidx[klane] += c * xr - rr3 * ukdx;
            shfidy[klane] += c * yr - rr3 * ukdy;
            shfidz[klane] += c * zr - rr3 * ukdz;


            c = rr5 * dot3(xr, yr, zr, uidx, uidy, uidz);
            fkdx += c * xr - rr3 * uidx;
            fkdy += c * yr - rr3 * uidy;
            fkdz += c * zr - rr3 * uidz;
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


} // generated by ComplexKernelBuilder (ck.py) 1.5.1


void sparse_precond_apply_cu2(const real (*rsd)[3], real (*zrsd)[3])
{
   const auto& st = *uspatial_v2_unit;
   const real off = switch_off(switch_usolve) + st.buffer;


   launch_k1s(nonblk, n, sparse_precond_cu3, //
              rsd, zrsd, polarity, n, udiag);

   int ngrid = get_grid_size(BLOCK_DIM);
   sparse_precond_cu4<<<ngrid, BLOCK_DIM, 0, nonblk>>>(
      st.n, TINKER_IMAGE_ARGS, off, st.si1.bit0, nwexclude, wexclude,
      wexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl, st.niak,
      st.iak, st.lst, rsd, zrsd, palpha, polarity);
}
}
