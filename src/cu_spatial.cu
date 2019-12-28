#include "box.h"
#include "launch.h"
#include "md.h"
#include "nblist.h"
#include "seq_image.h"
#include "seq_spatial_box.h"
#include "spatial.h"
#include "syntax/cu/ffsn.h"
#include "thrust_cache.h"
#include <thrust/extrema.h>
#include <thrust/remove.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/transform_reduce.h>
#include <thrust/transform_scan.h>


TINKER_NAMESPACE_BEGIN
namespace pbc {
#pragma acc routine seq
__device__
inline real3 c0_f0_triclinic(real xr, real yr, real zr, real3 ra, real3 rb,
                             real3 rc)
{
   real3 f;
   f.x = zr * ra.z + yr * ra.y + xr * ra.x;
   f.y = zr * rb.z + yr * rb.y;
   f.z = zr * rc.z;
   return f;
}


#pragma acc routine seq
__device__
inline real3 c0_f0_monoclinic(real xr, real yr, real zr, real3 ra, real3 rb,
                              real3 rc)
{
   real3 f;
   f.x = zr * ra.z + xr * ra.x;
   f.y = yr * rb.y;
   f.z = zr * rc.z;
   return f;
}


#pragma acc routine seq
__device__
inline real3 c0_f0_orthogonal(real xr, real yr, real zr, real3 ra, real3 rb,
                              real3 rc)
{
   real3 f;
   f.x = xr * ra.x;
   f.y = yr * rb.y;
   f.z = zr * rc.z;
   return f;
}


#pragma acc routine seq
__device__
inline real3 f0_f1(real3 f)
{
   f.x -= REAL_FLOOR(0.5f + f.x);
   f.y -= REAL_FLOOR(0.5f + f.y);
   f.z -= REAL_FLOOR(0.5f + f.z);
   return f;
}


#pragma acc routine seq
__device__
inline real3 f1_c1_triclinic(real3 f, real3 l1, real3 l2, real3 l3)
{
   f.x = f.z * l1.z + f.y * l1.y + f.x * l1.x;
   f.y = f.z * l2.z + f.y * l2.y;
   f.z = f.z * l3.z;
   return f;
}

#pragma acc routine seq
__device__
inline real3 f1_c1_monoclinic(real3 f, real3 l1, real3 l2, real3 l3)
{
   f.x = f.z * l1.z + f.x * l1.x;
   f.y = f.y * l2.y;
   f.z = f.z * l3.z;
   return f;
}

#pragma acc routine seq
__device__
inline real3 f1_c1_orthogonal(real3 f, real3 l1, real3 l2, real3 l3)
{
   f.x = f.x * l1.x;
   f.y = f.y * l2.y;
   f.z = f.z * l3.z;
   return f;
}
}

__device__
static inline real3 frac_general(real xr, real yr, real zr, real3 ra, real3 rb,
                                 real3 rc)
{
   using namespace pbc;
   if (ra.z == 0) {
      return f0_f1(c0_f0_orthogonal(xr, yr, zr, ra, rb, rc));
   } else if (ra.y == 0) {
      return f0_f1(c0_f0_monoclinic(xr, yr, zr, ra, rb, rc));
   } else {
      return f0_f1(c0_f0_triclinic(xr, yr, zr, ra, rb, rc));
   }
}


__device__
static inline real3 frac_image_general(real3 f, real3 l1, real3 l2, real3 l3)
{
   using namespace pbc;
   if (l1.z == 0) {
      return f1_c1_orthogonal(f, l1, l2, l3);
   } else if (l1.y == 0) {
      return f1_c1_monoclinic(f, l1, l2, l3);
   } else {
      return f1_c1_triclinic(f, l1, l2, l3);
   }
}


/**
 * \def frac
 * \ingroup macro
 * Calculate the fractional coordinates of (`xr, yr, zr`). The range of the
 * fractional coordinate is `[-1/2, 1/2)`.
 */
#ifndef frac
#   define frac(xr, yr, zr) frac_general(xr, yr, zr, recipa, recipb, recipc)
#endif
/**
 * \def frac_image
 * \ingroup macro
 */
#ifndef frac_image
#   define frac_image(f) frac_image_general(f, lvec1, lvec2, lvec3)
#endif


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
bool nearby_box0(int boxj, int px, int py, int pz, real3 lvec1, real3 lvec2,
                 real3 lvec3, real cutbuf2)
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
   // min(image(da), image(cb))
   real3 r = make_real3(0, 0, 0);
   if (2 <= ix && ix <= dimx - 2) {
      real da = ((real)ix + 1) / dimx;
      real cb = ((real)ix - 1) / dimx;
      da -= REAL_FLOOR(da + 0.5f);
      cb -= REAL_FLOOR(cb + 0.5f);
      r.x = REAL_MIN(REAL_ABS(da), REAL_ABS(cb));
   }
   if (2 <= iy && iy <= dimy - 2) {
      real da = ((real)iy + 1) / dimy;
      real cb = ((real)iy - 1) / dimy;
      da -= REAL_FLOOR(da + 0.5f);
      cb -= REAL_FLOOR(cb + 0.5f);
      r.y = REAL_MIN(REAL_ABS(da), REAL_ABS(cb));
   }
   if (2 <= iz && iz <= dimz - 2) {
      real da = ((real)iz + 1) / dimz;
      real cb = ((real)iz - 1) / dimz;
      da -= REAL_FLOOR(da + 0.5f);
      cb -= REAL_FLOOR(cb + 0.5f);
      r.z = REAL_MIN(REAL_ABS(da), REAL_ABS(cb));
   }
   r = frac_image(r);
   real r2 = r.x * r.x + r.y * r.y + r.z * r.z;
   return r2 <= cutbuf2;
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


extern "C" __global__
void spatial_bc(int n, int px, int py, int pz,
                Spatial::SortedAtom* restrict sorted, int* restrict boxnum,
                int* restrict nax, //
                const real* restrict x, const real* restrict y,
                const real* restrict z, TINKER_IMAGE_PARAMS, real cutbuf2,
                int ZERO_LBUF, real* restrict xold, real* restrict yold,
                real* restrict zold, //
                int nx, int* restrict nearby)
{
   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
        i += blockDim.x * gridDim.x) {
      real xr = x[i];
      real yr = y[i];
      real zr = z[i];
      if (!ZERO_LBUF) {
         xold[i] = xr;
         yold[i] = yr;
         zold[i] = zr;
      }
      real3 f = frac(xr, yr, zr);
      sorted[i].x = xr;       // B.2
      sorted[i].y = yr;       // B.2
      sorted[i].z = zr;       // B.2
      sorted[i].unsorted = i; // B.2
      int id = frac_to_box(px, py, pz, f.x, f.y, f.z);
      boxnum[i] = id;         // B.3
      atomicAdd(&nax[id], 1); // B.4
   }


   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < nx;
        i += blockDim.x * gridDim.x) {
      if (nearby_box0(i, px, py, pz, lvec1, lvec2, lvec3, cutbuf2))
         nearby[i] = i; // C.1 (close enough)
      else
         nearby[i] = -1; // C.1 (otherwise)
   }
}


extern "C" __global__
void spatial_e(int n, int nak, const int* restrict boxnum, int* xakf,
               const Spatial::SortedAtom* restrict sorted, TINKER_IMAGE_PARAMS)
{
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);
   const int prevlane = (ilane + WARP_SIZE - 1) & (WARP_SIZE - 1); // E.2


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


extern "C" __global__
void spatial_ghi(Spatial* restrict sp, int n, TINKER_IMAGE_PARAMS, real cutbuf2)
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
            real rik2 = imagen2(xr, yr, zr);
            int ilane_incl_j =
               (atomi < atomk && rik2 <= cutbuf2) ? 1 : 0; // I.5
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


extern "C" __global__
void spatial_update_sorted(int n, Spatial::SortedAtom* restrict sorted,
                           const real* restrict x, const real* restrict y,
                           const real* restrict z, TINKER_IMAGE_PARAMS)
{
   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
        i += blockDim.x * gridDim.x) {
      int ia = sorted[i].unsorted;
      real xr = x[ia];
      real yr = y[ia];
      real zr = z[ia];
      image(xr, yr, zr);
      sorted[i].x = xr;
      sorted[i].y = yr;
      sorted[i].z = zr;
   }
}
TINKER_NAMESPACE_END


TINKER_NAMESPACE_BEGIN
void spatial_data_update_sorted(SpatialUnit u)
{
   auto& st = *u;
   launch_k1s(nonblk, n, spatial_update_sorted, n, st.sorted, st.x, st.y, st.z,
              TINKER_IMAGE_ARGS);
}


void spatial_data_init_cu(SpatialUnit u)
{
   assert(u->rebuild == 1);
   u->rebuild = 0;
   const real cutbuf = u->cutoff + u->buffer;
   const real cutbuf2 = cutbuf * cutbuf;
   const real lbuf = u->buffer;
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
   auto policy = thrust::cuda::par(thrust_cache).on(nonblk);


   // B.1 D.1
   device_array::zero(false, nx + 1, ax_scan);
   // B.2 B.3 B.4 C.1
   const auto* lx = u->x;
   const auto* ly = u->y;
   const auto* lz = u->z;
   int ZERO_LBUF = (lbuf <= 0 ? 1 : 0);
   launch_k1s(nonblk, n, spatial_bc,                      //
              n, px, py, pz, sorted, boxnum, ax_scan + 1, //
              lx, ly, lz, TINKER_IMAGE_ARGS, cutbuf2, ZERO_LBUF, u->xold,
              u->yold, u->zold, //
              nx, nearby);
   // find max(nax) and compare to Spatial::BLOCK
   // ax_scan[0] == 0 can never be the maximum
   int level = 1 + floor_log2(nak - 1);
   int mnax;
   const int* mnaxptr = thrust::max_element(policy, ax_scan, ax_scan + 1 + nx);
   device_array::copyout(1, &mnax, mnaxptr, false);
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
      px = (level + 0) / 3;
      py = (level + 1) / 3;
      pz = (level + 2) / 3;
      nx = pow2(px + py + pz);
      nxk = (nx + Spatial::BLOCK - 1) / Spatial::BLOCK;

      device_array::allocate(nx, &nearby);
      device_array::allocate(nx + 1, &ax_scan);
      device_array::allocate(nak * nxk, &xkf);

      u.update_deviceptr(*u, false);

      device_array::zero(false, nx + 1, ax_scan);
      int ZERO_LBUF = (lbuf <= 0 ? 1 : 0);
      launch_k1s(nonblk, n, spatial_bc,                      //
                 n, px, py, pz, sorted, boxnum, ax_scan + 1, //
                 lx, ly, lz, TINKER_IMAGE_ARGS, cutbuf2, ZERO_LBUF, u->xold,
                 u->yold, u->zold, //
                 nx, nearby);
      mnaxptr = thrust::max_element(policy, ax_scan, ax_scan + 1 + nx);
      device_array::copyout(1, &mnax, mnaxptr, false);
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
   launch_k1s(nonblk, padded, spatial_e, n, nak, boxnum, xakf, sorted,
              TINKER_IMAGE_ARGS);
   // F.1
   xak_sum = thrust::transform_reduce(policy, xakf, xakf + nak, POPC(), 0,
                                      thrust::plus<int>());
   // F.2
   thrust::transform_exclusive_scan(policy, xakf, xakf + nak, xakf_scan, POPC(),
                                    0, thrust::plus<int>());
   if (xak_sum > xak_sum_cap) {
      device_array::deallocate(u->lst, u->iak);
      xak_sum_cap = xak_sum;
      size_t iak_size = near * xak_sum;            // F.3
      size_t lst_size = iak_size * Spatial::BLOCK; // F.4
      device_array::allocate(lst_size, &u->lst);
      device_array::allocate(iak_size, &u->iak);
   }
   // must update the device pointer to apply the changes in xak_sum
   u.update_deviceptr(*u, false);


   device_array::zero(false, near * xak_sum * Spatial::BLOCK,
                      u->lst);                // G.6
   device_array::zero(false, nak, naak);      // H.1
   device_array::zero(false, nak * nxk, xkf); // H.1
   launch_k1s(nonblk, padded, spatial_ghi, u.deviceptr(), n, TINKER_IMAGE_ARGS,
              cutbuf2);


   Int32* lst32 = (Int32*)u->lst;
   auto tup_begin =
      thrust::make_zip_iterator(thrust::make_tuple(u->iak, lst32));
   auto tup_end = thrust::make_zip_iterator(
      thrust::make_tuple(u->iak + near * xak_sum, lst32 + near * xak_sum));
   auto end2 = thrust::remove_if(policy, tup_begin, tup_end,
                                 IntInt32Pair::Int32IsZero());  // G.7
   u->niak = thrust::get<1>(end2.get_iterator_tuple()) - lst32; // G.7
   assert((thrust::get<0>(end2.get_iterator_tuple()) - u->iak) == u->niak);
   u.update_deviceptr(*u, false);
}
TINKER_NAMESPACE_END
