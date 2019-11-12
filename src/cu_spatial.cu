#include "box.h"
#include "cu_launch.h"
#include "mathfunc.cuh"
#include "md.h"
#include "nblist.h"
#include "seq_image.h"
#include "seq_spatial_box.h"
#include "spatial.h"
#include "thrust_cache.h"
#include <thrust/extrema.h>
#include <thrust/remove.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/transform_reduce.h>
#include <thrust/transform_scan.h>


TINKER_NAMESPACE_BEGIN
struct POPC
{
   __device__
   int operator()(int flag)
   {
      return __popc(flag);
   }
};


struct Int32
{
   long4 lx, ly, lz, lw;


   __device__
   static bool is_long4_zero(const long4& l)
   {
      return l.x == 0 && l.y == 0 && l.z == 0 && l.w == 0;
   }


   __device__
   static bool is_zero(const Int32& i32)
   {
      return is_long4_zero(i32.lx) && is_long4_zero(i32.ly) &&
         is_long4_zero(i32.lz) && is_long4_zero(i32.lw);
   }
};


struct IntInt32Pair
{
   struct Int32IsZero
   {
      __device__
      bool operator()(const thrust::tuple<int, Int32>& t)
      {
         return Int32::is_zero(thrust::get<1>(t));
      }
   };
};


__device__
bool nearby_box0(int boxj, int px, int py, int pz, const Box* restrict box,
                 real cutbuf)
{
   int dimx = 1 << px;
   int dimy = 1 << py;
   int dimz = 1 << pz;
   int ix, iy, iz;
   box_to_ixyz(ix, iy, iz, px, py, pz, boxj);

   // (a, b): (-0.5, a+1/dim)
   // (c, d): (a+ix/dim, c+1/dim)
   // da = a+(ix+1)/dim - a = (ix+1)/dim
   // cb = a+ix/dim - a-1/dim = (ix-1)/dim
   // min(imagen(da), imagen(cb))
   real rx = 0;
   real ry = 0;
   real rz = 0;
   if (2 <= ix && ix <= dimx - 2) {
      real da = ((real)ix + 1) / dimx;
      real cb = ((real)ix - 1) / dimx;
      da -= REAL_FLOOR(da + 0.5f);
      cb -= REAL_FLOOR(cb + 0.5f);
      rx = REAL_MIN(REAL_ABS(da), REAL_ABS(cb));
   }
   if (2 <= iy && iy <= dimy - 2) {
      real da = ((real)iy + 1) / dimy;
      real cb = ((real)iy - 1) / dimy;
      da -= REAL_FLOOR(da + 0.5f);
      cb -= REAL_FLOOR(cb + 0.5f);
      ry = REAL_MIN(REAL_ABS(da), REAL_ABS(cb));
   }
   if (2 <= iz && iz <= dimz - 2) {
      real da = ((real)iz + 1) / dimz;
      real cb = ((real)iz - 1) / dimz;
      da -= REAL_FLOOR(da + 0.5f);
      cb -= REAL_FLOOR(cb + 0.5f);
      rz = REAL_MIN(REAL_ABS(da), REAL_ABS(cb));
   }
   frac_image(rx, ry, rz, box);
   real r2 = rx * rx + ry * ry + rz * rz;
   if (r2 <= cutbuf * cutbuf)
      return true;
   else
      return false;
}


__device__
inline int offset_box(int nx, int ny, int nz, int ix1, int iy1, int iz1,
                      int offset)
{
   int dimx = (1 << nx);
   int dimy = (1 << ny);
   int dimz = (1 << nz);
   int ix, iy, iz;
   box_to_ixyz(ix, iy, iz, nx, ny, nz, offset);
   ix = (ix + ix1) & (dimx - 1);
   iy = (iy + iy1) & (dimy - 1);
   iz = (iz + iz1) & (dimz - 1);
   int id = ixyz_to_box(nx, ny, nz, ix, iy, iz);
   return id;
}


__global__
void spatial_bc(Spatial* restrict sp, const real* restrict x,
                const real* restrict y, const real* restrict z,
                const Box* restrict box, real cutbuf)
{
   int n = sp->n;
   int px = sp->px;
   int py = sp->py;
   int pz = sp->pz;
   auto* restrict sorted = sp->sorted;
   auto* restrict boxnum = sp->boxnum;
   auto* restrict nax = sp->ax_scan + 1; // D.2
   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
        i += blockDim.x * gridDim.x) {
      real xr = x[i];
      real yr = y[i];
      real zr = z[i];
      real fx, fy, fz;
      frac(fx, fy, fz, xr, yr, zr, box);
      sorted[i].x = xr;       // B.2
      sorted[i].y = yr;       // B.2
      sorted[i].z = zr;       // B.2
      sorted[i].unsorted = i; // B.2
      int id = frac_to_box(px, py, pz, fx, fy, fz);
      boxnum[i] = id;         // B.3
      atomicAdd(&nax[id], 1); // B.4
   }


   int nx = sp->nx;
   auto* restrict nearby = sp->nearby;
   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < nx;
        i += blockDim.x * gridDim.x) {
      if (nearby_box0(i, px, py, pz, box, cutbuf))
         nearby[i] = i; // C.1 (close enough)
      else
         nearby[i] = -1; // C.1 (otherwise)
   }
}


__global__
void spatial_e(Spatial* restrict sp)
{
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);
   const int prevlane = (ilane + WARP_SIZE - 1) & (WARP_SIZE - 1); // E.2


   const int n = sp->n;
   const int nak = sp->nak;
   const auto* restrict boxnum = sp->boxnum;
   auto* restrict xakf = sp->xakf;
   for (int iw = iwarp; iw < nak; iw += nwarp) {
      int atomi = iw * WARP_SIZE + ilane;
      int id1 = ((atomi < n) ? boxnum[atomi] : boxnum[n - 1]); // E.3
      int id0 = __shfl_sync(ALL_LANES, id1, prevlane);         // E.5
      int diff = (id0 == id1 ? 0 : 1);                         // E.1
      int flag = __ballot_sync(ALL_LANES, diff);               // E.6
      if (ilane == 0)
         xakf[iw] = (flag == 0 ? 1 : flag); // E.4
   }
}


__global__
void spatial_ghi(Spatial* restrict sp, const Box* restrict box, real cutbuf)
{
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);

   const int nak = sp->nak;
   const int px = sp->px;
   const int py = sp->py;
   const int pz = sp->pz;
   const int nxk = sp->nxk;
   const int near = sp->near;

   const auto* restrict boxnum = sp->boxnum;
   const auto* restrict xakf = sp->xakf;
   const auto* restrict xakf_scan = sp->xakf_scan;
   const auto* restrict nearby = sp->nearby;
   const auto* restrict begin = sp->ax_scan; // D.4
   const auto* restrict end = begin + 1;     // D.4

   auto* restrict iak = sp->iak;
   auto* restrict lst = sp->lst;
   auto* restrict naak = sp->naak;
   auto* restrict xkf = sp->xkf;

   for (int iw = iwarp; iw < nak; iw += nwarp) {
      int offset = xakf_scan[iw]; // F.5
      int flag = xakf[iw];        // E.7
      int nbox = __popc(flag);    // E.7

      auto* restrict iakbuf = iak + near * offset;             // G.4
      auto* restrict lstbuf = lst + near * offset * WARP_SIZE; // G.5
      auto* restrict ixkf = xkf + iw * nxk;                    // H.2
      const int atom_block_min = iw * WARP_SIZE;               // H.4
      for (int j = ilane; j < nbox * near; j += WARP_SIZE) {
         iakbuf[j] = iw;    // G.4
         int i0 = j / near; // the i-th least significant bit is i0 + 1
         int pos = ffsn(flag, i0 + 1) - 1;        // E.8
         int ibox = boxnum[iw * WARP_SIZE + pos]; // E.8
         int ix1, iy1, iz1;
         box_to_ixyz(ix1, iy1, iz1, px, py, pz, ibox);
         int j0 = nearby[j - i0 * near];
         int jbox = offset_box(px, py, pz, ix1, iy1, iz1, j0);
         // the (jbox%32)-th bit of the (jbox/32) flag will be set to 1
         int ii = jbox / WARP_SIZE;
         int jj = jbox & (WARP_SIZE - 1);
         int oldflag = atomicOr(&ixkf[ii], 1 << jj); // H.3
         // the atomicOr() will return the old value;
         // code in the following if body will only run
         // when the bit(ii,jj) gets set for the first time
         if ((oldflag & (1 << jj)) == 0) {
            // copy atoms in jbox to lstbuf
            int begin_i = begin[jbox];
            begin_i = max(atom_block_min + 1, begin_i);        // H.4
            int len = end[jbox] - begin_i;                     // H.5
            int start_pos = atomicAdd(&naak[iw], max(0, len)); // H.6
            // atomicAdd() will return the old value;
            // skip the loop if len is less than 1
            for (int kk = 0; kk < len; ++kk) {
               lstbuf[start_pos + kk] = begin_i + kk; // H.4
            }
         }
      }
   }


   real off2 = cutbuf * cutbuf;
   const int n = sp->n;
   const auto* restrict sorted = sp->sorted;
   for (int iw = iwarp; iw < nak; iw += nwarp) {
      int offset = xakf_scan[iw];
      const auto* restrict iakbuf = iak + near * offset;
      auto* restrict lstbuf = lst + near * offset * WARP_SIZE;
      int naak_coarse = naak[iw]; // I.1


      int start_pos = 0;
      int atomi;
      atomi = min(iakbuf[0] * WARP_SIZE + ilane, n - 1); // I.4
      real xi = sorted[atomi].x;
      real yi = sorted[atomi].y;
      real zi = sorted[atomi].z;
      int shatomk;
      real shx, shy, shz;
      int idx_max = (naak_coarse + WARP_SIZE - 1) / WARP_SIZE;
      idx_max *= WARP_SIZE; // I.2a
      for (int idx = ilane; idx < idx_max; idx += WARP_SIZE) {
         shatomk = lstbuf[idx]; // I.2b
         shx = sorted[shatomk].x;
         shy = sorted[shatomk].y;
         shz = sorted[shatomk].z;
         lstbuf[idx] = 0; // I.3


         int jflag = 0;
         for (int j = 0; j < WARP_SIZE; ++j) {
            int srclane = j;
            int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
            real xr = xi - __shfl_sync(ALL_LANES, shx, srclane);
            real yr = yi - __shfl_sync(ALL_LANES, shy, srclane);
            real zr = zi - __shfl_sync(ALL_LANES, shz, srclane);
            imagen(xr, yr, zr, box);
            real rik2 = xr * xr + yr * yr + zr * zr;
            int ilane_incl_j = (atomi < atomk && rik2 <= off2) ? 1 : 0; // I.5
            int incl_j = __ballot_sync(ALL_LANES, ilane_incl_j);
            if (incl_j)
               jflag |= (1 << j); // I.5
         }


         int njbit = __popc(jflag);
         int jth = ffsn(jflag, ilane + 1) - 1;
         int atomnb = __shfl_sync(ALL_LANES, shatomk, jth); // I.6a
         if (ilane < njbit)
            lstbuf[start_pos + ilane] = atomnb; // I.6b
         start_pos += njbit;
      }
   }
}
TINKER_NAMESPACE_END


TINKER_NAMESPACE_BEGIN
void spatial_data_init_cu(SpatialUnit u, NBListUnit nu)
{
   const real cutbuf = nu->cutoff + nu->buffer;
   const int& nak = u->nak;
   const int padded = nak * Spatial::BLOCK;
   int& px = u->px;
   int& py = u->py;
   int& pz = u->pz;
   int& nx = u->nx;
   int& nxk = u->nxk;
   int& near = u->near;
   int& xak_sum = u->xak_sum;
   int& xak_sum_cap = u->xak_sum_cap;
   int& niak = u->niak;


   auto*& sorted = u->sorted;
   auto*& boxnum = u->boxnum;
   auto*& naak = u->naak;
   auto*& xakf = u->xakf;
   auto*& xakf_scan = u->xakf_scan;
   auto*& nearby = u->nearby;
   auto*& ax_scan = u->ax_scan;
   auto*& xkf = u->xkf;


   // auto policy = thrust::device;
   auto policy = thrust::cuda::par(thrust_cache);


   // B.1 D.1
   device_array::zero(nx + 1, ax_scan);
   // B.2 B.3 B.4 C.1
   const auto* lx = nu->x;
   const auto* ly = nu->y;
   const auto* lz = nu->z;
   launch_kernel1(n, spatial_bc, u.deviceptr(), lx, ly, lz, box, cutbuf);
   // find max(nax) and compare to Spatial::BLOCK
   // ax_scan[0] == 0 can never be the maximum
   int level = 1 + floor_log2(nak - 1);
   int mnax;
   const int* mnaxptr = thrust::max_element(policy, ax_scan, ax_scan + 1 + nx);
   device_array::copyout(1, &mnax, mnaxptr);
   while (mnax > Spatial::BLOCK) {
      device_array::deallocate(nearby, ax_scan, xkf);

      int scale = (mnax - 1) / Spatial::BLOCK;
      // mnax / mnax-1 / scale / 2^p / p
      // 33   / 32     / 1     / 2   / 1
      // 64   / 63     / 1     / 2   / 1
      // 65   / 64     / 2     / 4   / 2
      // 128  / 127    / 3     / 4   / 2
      // 129  / 128    / 4     / 8   / 3
      int p = 1 + floor_log2(scale);
      level += p;
      px = (level + 2) / 3;
      py = (level + 1) / 3;
      pz = (level + 0) / 3;
      nx = ct::pow2(px + py + pz);
      nxk = (nx + Spatial::BLOCK - 1) / Spatial::BLOCK;

      device_array::allocate(nx, &nearby);
      device_array::allocate(nx + 1, &ax_scan);
      device_array::allocate(nak * nxk, &xkf);

      u.update_deviceptr(*u);

      device_array::zero(nx + 1, ax_scan);
      launch_kernel1(n, spatial_bc, u.deviceptr(), lx, ly, lz, box, cutbuf);
      mnaxptr = thrust::max_element(policy, ax_scan, ax_scan + 1 + nx);
      device_array::copyout(1, &mnax, mnaxptr);
   }
   // B.5
   thrust::stable_sort_by_key(policy, boxnum, boxnum + n, sorted);
   // C.2
   int* nearby_end = thrust::remove(policy, nearby, nearby + nx, -1);
   // C.3
   near = nearby_end - nearby;
   // D.2
   int* nax = ax_scan + 1;
   // D.3
   thrust::inclusive_scan(policy, nax, nax + nx, nax);


   // E
   launch_kernel1(padded, spatial_e, u.deviceptr());
   // F.1
   xak_sum = thrust::transform_reduce(policy, xakf, xakf + nak, POPC(), 0,
                                      thrust::plus<int>());
   // F.2
   thrust::transform_exclusive_scan(policy, xakf, xakf + nak, xakf_scan, POPC(),
                                    0, thrust::plus<int>());
   if (xak_sum > xak_sum_cap) {
      device_array::deallocate(u->iak, u->lst);
      xak_sum_cap = xak_sum;
      size_t iak_size = near * xak_sum;            // F.3
      size_t lst_size = iak_size * Spatial::BLOCK; // F.4
      // allocate iak and lst together
      device_array::allocate(iak_size + lst_size, &u->lst);
      u->iak = u->lst + lst_size;
   }
   // must update the device pointer to apply the changes in xak_sum
   u.update_deviceptr(*u);


   device_array::zero(near * xak_sum * Spatial::BLOCK, u->lst); // G.6
   device_array::zero(nak, naak);                               // H.1
   device_array::zero(nak * nxk, xkf);                          // H.1
   launch_kernel1(padded, spatial_ghi, u.deviceptr(), box, cutbuf);


   Int32* lst32 = (Int32*)u->lst;
   auto tup_begin =
      thrust::make_zip_iterator(thrust::make_tuple(u->iak, lst32));
   auto tup_end = thrust::make_zip_iterator(
      thrust::make_tuple(u->iak + near * xak_sum, lst32 + near * xak_sum));
   auto end2 = thrust::remove_if(policy, tup_begin, tup_end,
                                 IntInt32Pair::Int32IsZero());  // G.7
   u->niak = thrust::get<1>(end2.get_iterator_tuple()) - lst32; // G.7
   assert((thrust::get<0>(end2.get_iterator_tuple()) - u->iak) == u->niak);
   u.update_deviceptr(*u);
}
TINKER_NAMESPACE_END
