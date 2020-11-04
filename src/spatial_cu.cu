// Eventually thrust will drop c++11 support.
#define THRUST_IGNORE_DEPRECATED_CPP_DIALECT


#include "box.h"
#include "image.h"
#include "imagefc_cu.h"
#include "launch.h"
#include "md.h"
#include "nblist.h"
#include "spatial.h"
#include "syntax/cu/copysign.h"
#include "syntax/cu/ffsn.h"
#include "thrust_cache.h"
#include <thrust/extrema.h>
#include <thrust/remove.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/transform_reduce.h>
#include <thrust/transform_scan.h>


namespace tinker {
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


/**
 * \return
 * The original integer which has the smaller absolute value.
 */
__device__
inline int min_by_abs(int a, int b)
{
   return (abs(b) < abs(a)) ? b : a;
}


__device__
inline int ixyz_to_box(int ix, int iy, int iz, int px, int py, int pz)
{
   int id = (ix << (pz + py)) + (iy << pz) + iz;
   return id;
}


__device__
inline void box_to_ixyz(int& restrict ix, int& restrict iy, int& restrict iz,
                        int px, int py, int pz, int boxid)
{
   ix = boxid >> (pz + py);         // ix = boxid / (dimz*dimy)
   boxid &= ((1 << (pz + py)) - 1); // boxid = boxid % (dimz*dimy)
   iy = boxid >> pz;                // iy = boxid / dimz
   iz = boxid & ((1 << pz) - 1);    // iz = boxid % dimz
}


/**
 * \note
 * Fractional coordinates have to be taken care of by the image routine first.
 */
__device__
inline void frac_to_ixyz(int& restrict ix, int& restrict iy, int& restrict iz,
                         int px, int py, int pz, real fx, real fy, real fz)
{
   // cannot use iw = fw * (1 << pw) + (1 << pw) / 2;
   // with the implicit cast rules of C,
   // (int iw) <- (float fltw) <- fw * (1 << pw) + (1 << pw) / 2;
   // e.g., if pw = 3, fw = 0.49999998, fltw = 7.99999984, but 32-bit float
   // cannot represent this value accurately so that fltw <- 8.0f and will
   // give a wrong iw value 8 (should have been 7).
   // ix = fx * (1 << px) + (1 << px) / 2; // ix = (fx+half) * 2^px;
   // iy = fy * (1 << py) + (1 << py) / 2;
   // iz = fz * (1 << pz) + (1 << pz) / 2;

   // cannot use iw = fw * (1 << pw) + (1 << (pw - 1));
   // because pw may be 0, and 1 << -1 is undefined.

   ix = ((double)fx) * (1 << px) + (1 << px) / 2;
   iy = ((double)fy) * (1 << py) + (1 << py) / 2;
   iz = ((double)fz) * (1 << pz) + (1 << pz) / 2;
}


/**
 * \brief
 * Check the `ix, iy, iz` parameters of a given spatial box used for truncated
 * octahedron periodic boundaries. Update them with `ix', iy', iz'` of the box
 * image that is (at least partially) inside the truncated octahedron.
 *
 * The predicate is, if the fractional coordinate (from -0.5 to 0.5) of its
 * innear-most vertex is outside of the space defined by surfaces
 * `|x| + |y| + |z| = 3/4`, it is considered to be outside and needs updating.
 */
__device__
inline void ixyz_octahedron(int& restrict ix, int& restrict iy,
                            int& restrict iz, int px, int py, int pz)
{
   int qx = (1 << px);
   int qy = (1 << py);
   int qz = (1 << pz);
   int qx2 = qx / 2;
   int qy2 = qy / 2;
   int qz2 = qz / 2;


   // Translate by half box size so fractional coordinate can start from -0.5.
   int ix1 = ix - qx2;
   int iy1 = iy - qy2;
   int iz1 = iz - qz2;


   // The innear-most vertex.
   int ix2 = min_by_abs(ix1, ix1 + 1);
   int iy2 = min_by_abs(iy1, iy1 + 1);
   int iz2 = min_by_abs(iz1, iz1 + 1);
   // Fractional coordinate of the inner-most vertex; Range: [-0.5 to 0.5).
   real hx = (real)ix2 / qx;
   real hy = (real)iy2 / qy;
   real hz = (real)iz2 / qz;
   if (REAL_ABS(hx) + REAL_ABS(hy) + REAL_ABS(hz) > 0.75f) {
      // If iw2 == iw1+1, iw1-iw2 < 0, box is on (-0.5, 0)
      //    should move to right
      //    iw1 += qw2; iw1 -= -qw2; iw1 -= SIGN(qw2, iw1-iw2)
      // If iw2 == iw, iw1-iw2 == 0, box is on (0, 0.5)
      //    should move to left
      //    iw1 -= qw2; iw1 -= SIGN(qw2, 0); iw1 -= SIGN(qw2, iw1-iw2)
      //
      //    iw1 -= SIGN(qw2, iw1-iw2)
      ix1 -= INT_COPYSIGN(qx2, ix1 - ix2);
      iy1 -= INT_COPYSIGN(qy2, iy1 - iy2);
      iz1 -= INT_COPYSIGN(qz2, iz1 - iz2);
   }


   // Translate by half box size again.
   ix = ix1 + qx2;
   iy = iy1 + qy2;
   iz = iz1 + qz2;
}


__device__
inline bool nearby_box0(int px, int py, int pz, BoxShape box_shape, real3 lvec1,
                        real3 lvec2, real3 lvec3, int boxj, real cutbuf2)
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
   if (box_shape == OCT_BOX) {
      // I dunno how to calculate the shortest distance between two boxes with
      // truncated octahedron periodic boundary, so I just calculate all of the
      // 27 possible distances and pick the shortest one.
      // "waterglobe" with 19200 atoms is a good test for this PBC. Small
      // systems can only validate the image kernels, not spatial decomposition.
      r = make_real3(1, 1, 1);
      real r2 = 3;
      real rxs[3], rys[3], rzs[3];
      rxs[0] = (real)(ix - 1) / dimx;
      rxs[1] = (real)(ix) / dimx;
      rxs[2] = (real)(ix + 1) / dimx;
      rys[0] = (real)(iy - 1) / dimy;
      rys[1] = (real)(iy) / dimy;
      rys[2] = (real)(iy + 1) / dimy;
      rzs[0] = (real)(iz - 1) / dimz;
      rzs[1] = (real)(iz) / dimz;
      rzs[2] = (real)(iz + 1) / dimz;
      for (int i = 0; i < 3; ++i) {
         rxs[i] -= REAL_FLOOR(rxs[i] + 0.5f);
         rys[i] -= REAL_FLOOR(rys[i] + 0.5f);
         rzs[i] -= REAL_FLOOR(rzs[i] + 0.5f);
      }
      for (int xx = 0; xx < 3; ++xx) {
         for (int yy = 0; yy < 3; ++yy) {
            for (int zz = 0; zz < 3; ++zz) {
               real tx = rxs[xx];
               real ty = rys[yy];
               real tz = rzs[zz];
               if (REAL_ABS(tx) + REAL_ABS(ty) + REAL_ABS(tz) > 0.75f) {
                  tx -= REAL_SIGN(0.5f, tx);
                  ty -= REAL_SIGN(0.5f, ty);
                  tz -= REAL_SIGN(0.5f, tz);
               }
               real t2 = tx * tx + ty * ty + tz * tz;
               if (t2 < r2) {
                  r.x = tx;
                  r.y = ty;
                  r.z = tz;
                  r2 = t2;
               }
            }
         }
      }
   }
   r = ftoc(r);
   real r2 = r.x * r.x + r.y * r.y + r.z * r.z;
   return r2 <= cutbuf2;
}


__device__
inline int offset_box(int ix1, int iy1, int iz1, int px, int py, int pz,
                      int offset, BoxShape box_shape)
{
   int dimx = (1 << px);
   int dimy = (1 << py);
   int dimz = (1 << pz);
   int ix, iy, iz;
   box_to_ixyz(ix, iy, iz, px, py, pz, offset);
   ix = (ix + ix1) & (dimx - 1);
   iy = (iy + iy1) & (dimy - 1);
   iz = (iz + iz1) & (dimz - 1);
   if (box_shape == OCT_BOX)
      ixyz_octahedron(ix, iy, iz, px, py, pz);
   int id = ixyz_to_box(ix, iy, iz, px, py, pz);
   return id;
}


__global__
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
      sorted[i].x = xr;       // B.2
      sorted[i].y = yr;       // B.2
      sorted[i].z = zr;       // B.2
      sorted[i].unsorted = i; // B.2

      real3 f = imagectof(xr, yr, zr);

      // f.x, f.y, f.z coordinates have been "wrapped" into the PBC box,
      // but sorted[i] coordinates were not.

      int ix, iy, iz;
      frac_to_ixyz(ix, iy, iz, px, py, pz, f.x, f.y, f.z);
      ix &= ((1 << px) - 1);
      iy &= ((1 << py) - 1);
      iz &= ((1 << pz) - 1);
      if (box_shape == OCT_BOX)
         ixyz_octahedron(ix, iy, iz, px, py, pz);
      int id = ixyz_to_box(ix, iy, iz, px, py, pz);
      boxnum[i] = id;         // B.3
      atomicAdd(&nax[id], 1); // B.4
   }


   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < nx;
        i += blockDim.x * gridDim.x) {
      bool is_nearby0 =
         nearby_box0(px, py, pz, box_shape, lvec1, lvec2, lvec3, i, cutbuf2);
      if (is_nearby0)
         nearby[i] = i; // C.1 (close enough)
      else
         nearby[i] = -1; // C.1 (otherwise)
   }
}


__global__
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


__global__
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
         int jbox = offset_box(ix1, iy1, iz1, px, py, pz, j0, box_shape);
         // the (jbox%32)-th bit of the (jbox/32) flag will be set to 1
         int ii = jbox / WARP_SIZE;
         int jj = jbox & (WARP_SIZE - 1);
         int oldflag = atomicOr(&ixkf[ii], 1 << jj); // H.3
         // atomicOr() will return the old value;
         // code in the following if block will only run
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


__global__
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
}


namespace tinker {
void spatial_data_update_sorted(SpatialUnit u)
{
   auto& st = *u;
   launch_k1s(g::s0, n, spatial_update_sorted, n, st.sorted, st.x, st.y, st.z,
              TINKER_IMAGE_ARGS);
}


void spatial_data_init_cu(SpatialUnit u)
{
   u->fresh = -1; // 0xFFFFFFFF
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
   int& iak_cap = u->iak_cap;
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
   auto policy = thrust::cuda::par(thrust_cache).on(g::s0);


   // B.1 D.1
   darray::zero(g::q0, nx + 1, ax_scan);
   // B.2 B.3 B.4 C.1
   const auto* lx = u->x;
   const auto* ly = u->y;
   const auto* lz = u->z;
   int ZERO_LBUF = (lbuf <= 0 ? 1 : 0);
   launch_k1s(g::s0, n, spatial_bc,                       //
              n, px, py, pz, sorted, boxnum, ax_scan + 1, //
              lx, ly, lz, TINKER_IMAGE_ARGS, cutbuf2, ZERO_LBUF, u->xold,
              u->yold, u->zold, //
              nx, nearby);
   // find max(nax) and compare to Spatial::BLOCK
   // ax_scan[0] == 0 can never be the maximum
   int level = px + py + pz;
   int mnax;
   const int* mnaxptr = thrust::max_element(policy, ax_scan, ax_scan + 1 + nx);
   darray::copyout(g::q0, 1, &mnax, mnaxptr);
   wait_for(g::q0);
   while (mnax > Spatial::BLOCK) {
      darray::deallocate(nearby, ax_scan, xkf);

      int scale = (mnax - 1) / Spatial::BLOCK;
      // mnax / mnax-1 / scale / 2^p / p
      // 33   / 32     / 1     / 2   / 1
      // 64   / 63     / 1     / 2   / 1
      // 65   / 64     / 2     / 4   / 2
      // 128  / 127    / 3     / 4   / 2
      // 129  / 128    / 4     / 8   / 3
      int p = 1 + floor_log2(scale);
      level += p;
      spatial_cut(px, py, pz, level);
      nx = pow2(px + py + pz);
      nxk = (nx + Spatial::BLOCK - 1) / Spatial::BLOCK;

      darray::allocate(nx, &nearby);
      darray::allocate(nx + 1, &ax_scan);
      darray::allocate(nak * nxk, &xkf);

      u.update_deviceptr(*u, g::q0);

      darray::zero(g::q0, nx + 1, ax_scan);
      int ZERO_LBUF = (lbuf <= 0 ? 1 : 0);
      launch_k1s(g::s0, n, spatial_bc,                       //
                 n, px, py, pz, sorted, boxnum, ax_scan + 1, //
                 lx, ly, lz, TINKER_IMAGE_ARGS, cutbuf2, ZERO_LBUF, u->xold,
                 u->yold, u->zold, //
                 nx, nearby);
      mnaxptr = thrust::max_element(policy, ax_scan, ax_scan + 1 + nx);
      darray::copyout(g::q0, 1, &mnax, mnaxptr);
      wait_for(g::q0);
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
   launch_k1s(g::s0, padded, spatial_e, n, nak, boxnum, xakf, sorted,
              TINKER_IMAGE_ARGS);
   // F.1
   xak_sum = thrust::transform_reduce(policy, xakf, xakf + nak, POPC(), 0,
                                      thrust::plus<int>());
   // F.2
   thrust::transform_exclusive_scan(policy, xakf, xakf + nak, xakf_scan, POPC(),
                                    0, thrust::plus<int>());
   size_t iak_size = near * xak_sum;            // F.3
   size_t lst_size = iak_size * Spatial::BLOCK; // F.4
   if (iak_size > iak_cap) {
      iak_cap = iak_size;
      darray::deallocate(u->lst, u->iak);
      darray::allocate(lst_size, &u->lst);
      darray::allocate(iak_size, &u->iak);
   }
   // must update the device pointer to apply the changes in xak_sum
   u.update_deviceptr(*u, g::q0);


   darray::zero(g::q0, near * xak_sum * Spatial::BLOCK, u->lst); // G.6
   darray::zero(g::q0, nak, naak);                               // H.1
   darray::zero(g::q0, nak * nxk, xkf);                          // H.1
   launch_k1s(g::s0, padded, spatial_ghi, u.deviceptr(), n, TINKER_IMAGE_ARGS,
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
   u.update_deviceptr(*u, g::q0);
}
}


#if 0
namespace tinker {
namespace spatial_v2 {
__device__
inline void box_to_ixyz(int& restrict ix, int& restrict iy, int& restrict iz,
                        int px, int py, int pz, int boxid)
{
   ix = 0;
   iy = 0;
   iz = 0;
   for (int i = 0; i < pz; ++i) {
      int zid = (boxid >> (2 * i + 0)) & (1 << i);
      int yid = (boxid >> (2 * i + 1)) & (1 << i);
      int xid = (boxid >> (2 * i + 2)) & (1 << i);
      iz += zid;
      iy += yid;
      ix += xid;
   }
   for (int i = pz; i < py; ++i) {
      int yid = (boxid >> (2 * i + 0)) & (1 << i);
      int xid = (boxid >> (2 * i + 1)) & (1 << i);
      iy += yid;
      ix += xid;
   }
   for (int i = py; i < px; ++i) {
      int xid = (boxid >> (2 * i + 0)) & (1 << i);
      ix += xid;
   }
}


__device__
inline int ixyz_to_box(int ix, int iy, int iz, int px, int py, int pz)
{
   int id = 0;
   for (int i = 0; i < pz; ++i) {
      int zid = (iz & (1 << i)) << (2 * i + 0);
      int yid = (iy & (1 << i)) << (2 * i + 1);
      int xid = (ix & (1 << i)) << (2 * i + 2);
      id += (zid + yid + xid);
   }
   for (int i = pz; i < py; ++i) {
      int yid = (iy & (1 << i)) << (2 * i + 0);
      int xid = (ix & (1 << i)) << (2 * i + 1);
      id += (yid + xid);
   }
   for (int i = py; i < px; ++i) {
      int xid = (ix & (1 << i)) << (2 * i + 0);
      id += xid;
   }
   return id;
}
}
}
#endif


//====================================================================//


#include "seq_triangle.h"
#include "spatial2.h"


namespace tinker {
namespace {
using coord_t = int;


template <size_t n>
__device__
inline void AxesToTranspose(coord_t (&x)[n], int b)
{
   coord_t M = 1 << (b - 1), t;
   for (coord_t Q = M; Q > 1; Q >>= 1) {
      coord_t P = Q - 1;
      #pragma unroll
      for (int i = 0; i < n; ++i) {
         if (x[i] & Q) {
            x[0] ^= P; // invert
         } else {
            t = (x[0] ^ x[i]) & P;
            x[0] ^= t;
            x[i] ^= t; // exchange
         }
      }
   }


   #pragma unroll
   for (int i = 1; i < n; ++i)
      x[i] ^= x[i - 1];


   t = 0;
   for (coord_t Q = M; Q > 1; Q >>= 1)
      if (x[n - 1] & Q)
         t ^= Q - 1;


   #pragma unroll
   for (int i = 0; i < n; ++i)
      x[i] ^= t;
}


template <size_t n>
__device__
inline int TransposeToIndex(coord_t (&x)[n], int b)
{
   int val = 0;
   for (int ib = b - 1; ib >= 0; --ib) {
      #pragma unroll
      for (int in = 0; in < n; ++in) {
         val <<= 1;
         val += (x[in] >> ib & 1);
      }
   }
   return val;
}
}


__global__
void spatial2_step1(int n, int pz, int2* restrict b2num, //
                    const real* restrict x, const real* restrict y,
                    const real* restrict z, TINKER_IMAGE_PARAMS, int nakpk,
                    int* restrict akpf)
{
   // i = unsorted atom number
   // b2num[i] = [box number][unsorted atom number]
   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
        i += blockDim.x * gridDim.x) {
      real xr = x[i];
      real yr = y[i];
      real zr = z[i];
      real3 f = imagectof(xr, yr, zr);


      int ix, iy, iz;
      frac_to_ixyz(ix, iy, iz, pz, pz, pz, f.x, f.y, f.z);
      // Due to the limited precision, f.z (and f.x, f.y) may turn out to be
      // +0.5, which should have been in the range of [-0.5, +0.5). Thus iz may
      // be (2**pz)/2, which is also out of range. However, the Hilbert curve
      // algorithm used here works so well that we don't even need to call
      // "image" for the fractional coordinates. Here we still "filter" iz as
      // follows because we want it to be "in-range" for the truncated
      // octahedron box.
      int px = pz, py = pz;
      ix &= ((1 << px) - 1);
      iy &= ((1 << py) - 1);
      iz &= ((1 << pz) - 1);
      if (box_shape == OCT_BOX)
         ixyz_octahedron(ix, iy, iz, pz, pz, pz);
      coord_t ixyz[3] = {ix, iy, iz};
      AxesToTranspose(ixyz, pz);
      int id = TransposeToIndex(ixyz, pz);
      b2num[i] = make_int2(id, i); // B.1
      // For debugging purpose, uncomment the next line to disable the sorting
      // in the next step, so that sorted[i].unsorted == i.
      // b2num[i] = make_int2(i, i);
   }


   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < nakpk;
        i += blockDim.x * gridDim.x) {
      akpf[i] = 0; // B.2
   }
}


__global__
void spatial2_update_sorted(int n, Spatial::SortedAtom* restrict sorted,
                            const real* restrict x, const real* restrict y,
                            const real* restrict z, //
                            TINKER_IMAGE_PARAMS, real cut,
                            Spatial2::Center* restrict akc,
                            Spatial2::Center* restrict half)
{
   real xbox, ybox, zbox;
   xbox = lvec1.x * lvec1.x + lvec2.x * lvec2.x + lvec3.x * lvec3.x;
   ybox = lvec1.y * lvec1.y + lvec2.y * lvec2.y + lvec3.y * lvec3.y;
   zbox = lvec1.z * lvec1.z + lvec2.z * lvec2.z + lvec3.z * lvec3.z;
   xbox = REAL_SQRT(xbox);
   ybox = REAL_SQRT(ybox);
   zbox = REAL_SQRT(zbox);


   real xr, yr, zr, r, r2;


   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
        i += blockDim.x * gridDim.x) {
      int atomi = sorted[i].unsorted;
      xr = x[atomi];
      yr = y[atomi];
      zr = z[atomi];
      sorted[i].x = xr;
      sorted[i].y = yr;
      sorted[i].z = zr;


      // C.4 mid point
      int tmask = __activemask();
      real x0 = __shfl_sync(tmask, xr, 0);
      real y0 = __shfl_sync(tmask, yr, 0);
      real z0 = __shfl_sync(tmask, zr, 0);
      // vector r(i)-r(0)
      xr -= x0;
      yr -= y0;
      zr -= z0;
      image2(xr, yr, zr);
      real xmin = xr;
      real ymin = yr;
      real zmin = zr;
      real xmax = xr;
      real ymax = yr;
      real zmax = zr;
      #pragma unroll
      for (int ii = 1; ii < WARP_SIZE; ii *= 2) {
         xmin = REAL_MIN(xmin, __shfl_xor_sync(tmask, xmin, ii));
         ymin = REAL_MIN(ymin, __shfl_xor_sync(tmask, ymin, ii));
         zmin = REAL_MIN(zmin, __shfl_xor_sync(tmask, zmin, ii));
         xmax = REAL_MAX(xmax, __shfl_xor_sync(tmask, xmax, ii));
         ymax = REAL_MAX(ymax, __shfl_xor_sync(tmask, ymax, ii));
         zmax = REAL_MAX(zmax, __shfl_xor_sync(tmask, zmax, ii));
      }
      real xc = (xmin + xmax) / 2;
      real yc = (ymin + ymax) / 2;
      real zc = (zmin + zmax) / 2;
      // C.4 radius
      r2 = (xr - xc) * (xr - xc);
      r2 += (yr - yc) * (yr - yc);
      r2 += (zr - zc) * (zr - zc);
      r = REAL_SQRT(r2);
      #pragma unroll
      for (int ii = 1; ii < WARP_SIZE; ii *= 2) {
         r = REAL_MAX(r, __shfl_xor_sync(tmask, r, ii));
      }
      // C.4 half size
      real xh = (xmax - xmin) / 2;
      real yh = (ymax - ymin) / 2;
      real zh = (zmax - zmin) / 2;
      // C.4 local flag
      // Different from list construction: instead of `cutoff+buffer`,
      // only `cutoff` is needed.
      // if i-block has a small distribution in space:
      // xbox / 2 - xhalf > cut ==> xbox > 2*(xhalf+cut)
      bool ilocal = xbox > 2 * (xh + cut) //
         and ybox > 2 * (yh + cut)        //
         and zbox > 2 * (zh + cut);


      int iblock = i / WARP_SIZE;
      int ilane = threadIdx.x & (WARP_SIZE - 1);
      if (ilane == 0) {
         akc[iblock].x = xc + x0;
         akc[iblock].y = yc + y0;
         akc[iblock].z = zc + z0;
         akc[iblock].w = ilocal;
         half[iblock].x = xh;
         half[iblock].y = yh;
         half[iblock].z = zh;
         half[iblock].w = r;
      }
   }
}


__global__
void spatial2_step2(int n, Spatial::SortedAtom* restrict sorted,
                    int* restrict bnum, int2* restrict b2num,
                    const real* restrict x, const real* restrict y,
                    const real* restrict z, int ZERO_LBUF, real* restrict xold,
                    real* restrict yold, real* restrict zold, //
                    TINKER_IMAGE_PARAMS, real cutbuf,
                    Spatial2::Center* restrict akc,
                    Spatial2::Center* restrict half)
{
   real xbox, ybox, zbox;
   xbox = lvec1.x * lvec1.x + lvec2.x * lvec2.x + lvec3.x * lvec3.x;
   ybox = lvec1.y * lvec1.y + lvec2.y * lvec2.y + lvec3.y * lvec3.y;
   zbox = lvec1.z * lvec1.z + lvec2.z * lvec2.z + lvec3.z * lvec3.z;
   xbox = REAL_SQRT(xbox);
   ybox = REAL_SQRT(ybox);
   zbox = REAL_SQRT(zbox);


   real xr, yr, zr, r, r2;


   // b2num has been sorted by the box number
   // i: sorted index
   // b2num[i] = [box number][unsorted atom number]
   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
        i += blockDim.x * gridDim.x) {
      int atomi = b2num[i].y;
      xr = x[atomi];
      yr = y[atomi];
      zr = z[atomi];
      if (!ZERO_LBUF) {
         xold[atomi] = xr;
         yold[atomi] = yr;
         zold[atomi] = zr;
      }
      sorted[i].x = xr;
      sorted[i].y = yr;
      sorted[i].z = zr;
      sorted[i].unsorted = atomi;
      bnum[atomi] = i;            // C.1
      b2num[i] = make_int2(0, 0); // C.3


      // C.4 mid point
      int tmask = __activemask();
      real x0 = __shfl_sync(tmask, xr, 0);
      real y0 = __shfl_sync(tmask, yr, 0);
      real z0 = __shfl_sync(tmask, zr, 0);
      // vector r(i)-r(0)
      xr -= x0;
      yr -= y0;
      zr -= z0;
      image2(xr, yr, zr);
      real xmin = xr;
      real ymin = yr;
      real zmin = zr;
      real xmax = xr;
      real ymax = yr;
      real zmax = zr;
      #pragma unroll
      for (int ii = 1; ii < WARP_SIZE; ii *= 2) {
         xmin = REAL_MIN(xmin, __shfl_xor_sync(tmask, xmin, ii));
         ymin = REAL_MIN(ymin, __shfl_xor_sync(tmask, ymin, ii));
         zmin = REAL_MIN(zmin, __shfl_xor_sync(tmask, zmin, ii));
         xmax = REAL_MAX(xmax, __shfl_xor_sync(tmask, xmax, ii));
         ymax = REAL_MAX(ymax, __shfl_xor_sync(tmask, ymax, ii));
         zmax = REAL_MAX(zmax, __shfl_xor_sync(tmask, zmax, ii));
      }
      real xc = (xmin + xmax) / 2;
      real yc = (ymin + ymax) / 2;
      real zc = (zmin + zmax) / 2;
      // C.4 radius
      r2 = (xr - xc) * (xr - xc);
      r2 += (yr - yc) * (yr - yc);
      r2 += (zr - zc) * (zr - zc);
      r = REAL_SQRT(r2);
      #pragma unroll
      for (int ii = 1; ii < WARP_SIZE; ii *= 2) {
         r = REAL_MAX(r, __shfl_xor_sync(tmask, r, ii));
      }
      // C.4 half size
      real xh = (xmax - xmin) / 2;
      real yh = (ymax - ymin) / 2;
      real zh = (zmax - zmin) / 2;
      // C.4 local flag
      // if i-block has a small distribution in space:
      // xbox / 2 - xhalf > cutbuf ==> xbox > 2*(xhalf+cutbuf)
      bool ilocal = xbox > 2 * (xh + cutbuf) //
         and ybox > 2 * (yh + cutbuf)        //
         and zbox > 2 * (zh + cutbuf);


      int iblock = i / WARP_SIZE;
      int ilane = threadIdx.x & (WARP_SIZE - 1);
      if (ilane == 0) {
         akc[iblock].x = xc + x0;
         akc[iblock].y = yc + y0;
         akc[iblock].z = zc + z0;
         akc[iblock].w = ilocal;
         half[iblock].x = xh;
         half[iblock].y = yh;
         half[iblock].z = zh;
         half[iblock].w = r;
      }
   }
}


__device__
inline void spatial2_step3_atomicOr(int x0, int y0, int* akpf, int* sum_nakpl)
{
   int x = max(x0, y0);
   int y = min(x0, y0);
   int f = xy_to_tri(x, y);
   int j = f / WARP_SIZE;
   int k = f & (WARP_SIZE - 1); // f % 32
   int mask = 1 << k;
   int oldflag = atomicOr(&akpf[j], mask);
   int oldkbit = oldflag & mask;
   if (oldkbit == 0) {
      atomicAdd(sum_nakpl, 1);
   }
}


__global__
void spatial2_step3(int nak, int* restrict akpf, int* nakpl_ptr0, //
                    const int* restrict bnum, int nstype,         //
                    int ns1, int (*restrict js1)[2],              //
                    int ns2, int (*restrict js2)[2],              //
                    int ns3, int (*restrict js3)[2],              //
                    int ns4, int (*restrict js4)[2])
{
   // D.1 Pairwise flag for (block i - block i) is always set.
   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < nak;
        i += blockDim.x * gridDim.x) {
      spatial2_step3_atomicOr(i, i, akpf, nakpl_ptr0);
   }


   // pairwise flag
   int maxns = -1;
   maxns = max(maxns, ns1);
   maxns = max(maxns, ns2);
   maxns = max(maxns, ns3);
   maxns = max(maxns, ns4);
   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < maxns;
        i += blockDim.x * gridDim.x) {
      int x0, y0;
      if (nstype >= 1 and i < ns1) {
         x0 = bnum[js1[i][0]] / WARP_SIZE;
         y0 = bnum[js1[i][1]] / WARP_SIZE;
         spatial2_step3_atomicOr(x0, y0, akpf, nakpl_ptr0);
      }
      if (nstype >= 2 and i < ns2) {
         x0 = bnum[js2[i][0]] / WARP_SIZE;
         y0 = bnum[js2[i][1]] / WARP_SIZE;
         spatial2_step3_atomicOr(x0, y0, akpf, nakpl_ptr0);
      }
      if (nstype >= 3 and i < ns3) {
         x0 = bnum[js3[i][0]] / WARP_SIZE;
         y0 = bnum[js3[i][1]] / WARP_SIZE;
         spatial2_step3_atomicOr(x0, y0, akpf, nakpl_ptr0);
      }
      if (nstype >= 4 and i < ns4) {
         x0 = bnum[js4[i][0]] / WARP_SIZE;
         y0 = bnum[js4[i][1]] / WARP_SIZE;
         spatial2_step3_atomicOr(x0, y0, akpf, nakpl_ptr0);
      }
   }
}


__global__
void spatial2_step4(int nakpk, int* restrict nakpl_ptr1,
                    const int* restrict akpf, int* restrict iakpl,
                    int* restrict iakpl_rev,   //
                    int cap_nakpl, int nstype, //
                    unsigned int* restrict s1b0, unsigned int* restrict s2b0,
                    unsigned int* restrict s3b0, unsigned int* restrict s4b0)
{
   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < nakpk;
        i += blockDim.x * gridDim.x) {
      int flag = akpf[i];
      int count = __popc(flag);
      int base = atomicAdd(nakpl_ptr1, count);
      int c = 0;
      while (flag) {
         int j = __ffs(flag) - 1;
         flag &= (flag - 1);
         int tri = WARP_SIZE * i + j;
         iakpl[base + c] = tri;
         iakpl_rev[tri] = base + c;
         ++c;
      }
   }


   // zero out bit0 and bit1
   for (int i = threadIdx.x + blockIdx.x * blockDim.x;
        i < WARP_SIZE * cap_nakpl; i += blockDim.x * gridDim.x) {
      if (nstype >= 1) {
         s1b0[i] = 0;
      }
      if (nstype >= 2) {
         s2b0[i] = 0;
      }
      if (nstype >= 3) {
         s3b0[i] = 0;
      }
      if (nstype >= 4) {
         s4b0[i] = 0;
      }
   }
}


__device__
void spatial2_step5_bits(int x0, int y0, unsigned int* bit0,
                         const int* iakpl_rev)
{
   int x, y;
   int bx, by, ax, ay;
   int f, fshort, pos;
   x = max(x0, y0);
   y = min(x0, y0);
   bx = x / WARP_SIZE;
   ax = x & (WARP_SIZE - 1);
   by = y / WARP_SIZE;
   ay = y & (WARP_SIZE - 1);
   f = xy_to_tri(bx, by);
   fshort = iakpl_rev[f];
   pos = WARP_SIZE * fshort + ax;
   atomicOr(&bit0[pos], 1 << ay);
   if (bx == by) {
      pos = WARP_SIZE * fshort + ay;
      atomicOr(&bit0[pos], 1 << ax);
   }
}


template <class IMG>
__global__
void spatial2_step5(const int* restrict bnum, const int* iakpl_rev, int nstype,
                    Spatial2::ScaleInfo si1, Spatial2::ScaleInfo si2,
                    Spatial2::ScaleInfo si3, Spatial2::ScaleInfo si4, //
                    int* restrict dev_niak, int* restrict iak,
                    int* restrict lst,                                //
                    int n, int nak, real cutbuf, TINKER_IMAGE_PARAMS, //
                    const int* restrict akpf,
                    const Spatial::SortedAtom* restrict sorted,
                    const Spatial2::Center* restrict akc,
                    const Spatial2::Center* restrict half)
{
   int maxns = -1;
   maxns = max(maxns, si1.ns);
   maxns = max(maxns, si2.ns);
   maxns = max(maxns, si3.ns);
   maxns = max(maxns, si4.ns);
   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < maxns;
        i += blockDim.x * gridDim.x) {
      int x0, y0;
      if (nstype >= 1 and i < si1.ns) {
         auto& si = si1;
         x0 = bnum[si.js[i][0]];
         y0 = bnum[si.js[i][1]];
         spatial2_step5_bits(x0, y0, si.bit0, iakpl_rev);
      }
      if (nstype >= 2 and i < si2.ns) {
         auto& si = si2;
         x0 = bnum[si.js[i][0]];
         y0 = bnum[si.js[i][1]];
         spatial2_step5_bits(x0, y0, si.bit0, iakpl_rev);
      }
      if (nstype >= 3 and i < si3.ns) {
         auto& si = si3;
         x0 = bnum[si.js[i][0]];
         y0 = bnum[si.js[i][1]];
         spatial2_step5_bits(x0, y0, si.bit0, iakpl_rev);
      }
      if (nstype >= 4 and i < si4.ns) {
         auto& si = si4;
         x0 = bnum[si.js[i][0]];
         y0 = bnum[si.js[i][1]];
         spatial2_step5_bits(x0, y0, si.bit0, iakpl_rev);
      }
   }


   int ithread, iwarp, nwarp, ilane;
   ithread = threadIdx.x + blockIdx.x * blockDim.x;
   iwarp = ithread / WARP_SIZE;
   nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   ilane = threadIdx.x & (WARP_SIZE - 1);


   constexpr int LEN_LSTBUF = BLOCK_DIM;
   __shared__ int lstbuf[LEN_LSTBUF * (BLOCK_DIM / WARP_SIZE)];
   int* buffer = &lstbuf[LEN_LSTBUF * (threadIdx.x / WARP_SIZE)];
   int lanemask = (1 << ilane) - 1;


   real xr, yr, zr, r2;
   real rlimit, rlimit2;
   real cutbuf2 = cutbuf * cutbuf;


   // Every warp loads a block of atoms as "i-block", denoted by "wy".
   for (int wy = iwarp; wy < nak - 1; wy += nwarp) {
      int atomi = wy * WARP_SIZE + ilane;
      real xi = sorted[atomi].x;
      real yi = sorted[atomi].y;
      real zi = sorted[atomi].z;
      real hxi = half[wy].x;
      real hyi = half[wy].y;
      real hzi = half[wy].z;
      real rad2 = half[wy].w;
      real cxi = akc[wy].x;
      real cyi = akc[wy].y;
      real czi = akc[wy].z;
      bool ilocal = akc[wy].w != 0;
      if (ilocal) {
         xr = xi - cxi;
         yr = yi - cyi;
         zr = zi - czi;
         IMG::img2(xr, yr, zr, TINKER_IMAGE_LVEC_ARGS, TINKER_IMAGE_RECIP_ARGS);
         xi = xr + cxi;
         yi = yr + cyi;
         zi = zr + czi;
      }
      // Number of k-neighbors found.
      int nknb = 0;


      // Every thread in the warp loads the center and size a bounding box,
      // denoted by "wx".
      for (int wx0 = wy + 1; wx0 < nak; wx0 += WARP_SIZE) {
         int wx = wx0 + ilane;
         bool calcwx = wx < nak; // wx cannot exceed nak-1.


         // If this block pair was recorded, we will skip it.
         if (calcwx) {
            int iw, iwa, iwb;
            iw = xy_to_tri(wx, wy);
            iwa = iw / WARP_SIZE;
            iwb = iw & (WARP_SIZE - 1);
            if (akpf[iwa] & (1 << iwb))
               calcwx = false;
         }


         // If this block pair was not recorded, but the blocks are far apart,
         // we will skip it.
         if (calcwx) {
            real cxk = akc[wx].x;
            real cyk = akc[wx].y;
            real czk = akc[wx].z;
            real rad1 = half[wx].w;
            rlimit = cutbuf + rad1 + rad2;
            rlimit2 = rlimit * rlimit;
            xr = cxk - cxi;
            yr = cyk - cyi;
            zr = czk - czi;
            r2 = IMG::img2(xr, yr, zr, TINKER_IMAGE_LVEC_ARGS,
                           TINKER_IMAGE_RECIP_ARGS);
            if (r2 > rlimit2)
               calcwx = false;
         }


         // Check if the bounding boxes of two blocks overlap.
         if (calcwx) {
            rlimit2 = cutbuf2;
            xr = REAL_ABS(xr) - (half[wx].x + hxi);
            yr = REAL_ABS(yr) - (half[wx].y + hyi);
            zr = REAL_ABS(zr) - (half[wx].z + hzi);
            xr = REAL_MAX((real)0, xr);
            yr = REAL_MAX((real)0, yr);
            zr = REAL_MAX((real)0, zr);
            r2 = xr * xr + yr * yr + zr * zr;
            if (r2 > rlimit2)
               calcwx = false;
         }


         // Only loop over the k-blocks that are in-use.
         int wxflag = __ballot_sync(ALL_LANES, calcwx);
         while (wxflag) {
            // Get the rightmost 1 digit.
            int jlane = __ffs(wxflag) - 1;
            // Clear the rightmost 1 digit.
            wxflag &= (wxflag - 1);


            wx = wx0 + jlane;
            real cxk = akc[wx].x;
            real cyk = akc[wx].y;
            real czk = akc[wx].z;
            real rad1 = half[wx].w;
            int k = wx * WARP_SIZE + ilane;
            int atomk = min(k, n - 1);
            real xk = sorted[atomk].x;
            real yk = sorted[atomk].y;
            real zk = sorted[atomk].z;
            if (ilocal) {
               xr = xk - cxi;
               yr = yk - cyi;
               zr = zk - czi;
               IMG::img2(xr, yr, zr, TINKER_IMAGE_LVEC_ARGS,
                         TINKER_IMAGE_RECIP_ARGS);
               xk = xr + cxi;
               yk = yr + cyi;
               zk = zr + czi;
            }


            // Check i-atoms and k-center.
            rlimit = rad1 + cutbuf;
            rlimit2 = rlimit * rlimit;
            xr = cxk - xi;
            yr = cyk - yi;
            zr = czk - zi;
            r2 = IMG::img2(xr, yr, zr, TINKER_IMAGE_LVEC_ARGS,
                           TINKER_IMAGE_RECIP_ARGS);
            // Only inlcude the "i-atoms" that are close to k-center.
            int iflag = __ballot_sync(ALL_LANES, r2 <= rlimit2);


            int includek = 0;
            if (iflag != 0) {
               if (ilocal) {
                  while (iflag) {
                     int j = __ffs(iflag) - 1;
                     iflag &= (iflag - 1);
                     xr = __shfl_sync(ALL_LANES, xi, j) - xk;
                     yr = __shfl_sync(ALL_LANES, yi, j) - yk;
                     zr = __shfl_sync(ALL_LANES, zi, j) - zk;
                     r2 = xr * xr + yr * yr + zr * zr;
                     if (r2 <= cutbuf2)
                        includek = 1;
                  }
               } else {
                  while (iflag) {
                     int j = __ffs(iflag) - 1;
                     iflag &= (iflag - 1);
                     xr = __shfl_sync(ALL_LANES, xi, j) - xk;
                     yr = __shfl_sync(ALL_LANES, yi, j) - yk;
                     zr = __shfl_sync(ALL_LANES, zi, j) - zk;
                     r2 = IMG::img2(xr, yr, zr, TINKER_IMAGE_LVEC_ARGS,
                                    TINKER_IMAGE_RECIP_ARGS);
                     if (r2 <= cutbuf2)
                        includek = 1;
                  }
               }
            }
            includek = includek and (k < n);
            int kflag = __ballot_sync(ALL_LANES, includek);


            if (includek) {
               int pos = nknb + __popc(kflag & lanemask);
               buffer[pos] = k;
            }
            nknb += __popc(kflag);
            if (nknb > LEN_LSTBUF - WARP_SIZE) {
               // Buffer is almost full, and there may not be enough space for
               // the next warp of k-atoms.
               int incr = nknb / WARP_SIZE;
               int pos;
               if (incr > 0) {
                  if (ilane == 0)
                     pos = atomicAdd(dev_niak, incr);
                  pos = __shfl_sync(ALL_LANES, pos, 0);
                  if (pos + incr <= nak * Spatial2::LSTCAP) {
                     if (ilane < incr)
                        iak[pos + ilane] = wy;
                     for (int i = 0; i < incr; ++i) {
                        int bufp = i * WARP_SIZE + ilane;
                        lst[(pos + i) * WARP_SIZE + ilane] = buffer[bufp];
                     }
                  }
                  if (incr < LEN_LSTBUF / WARP_SIZE)
                     buffer[ilane] = buffer[incr * WARP_SIZE + ilane];
                  nknb -= incr * WARP_SIZE;
               }
            }
         } // end while (wxflag)
      }


      if (nknb > 0) {
         int incr = (nknb + WARP_SIZE - 1) / WARP_SIZE;
         int pos;
         if (ilane == 0)
            pos = atomicAdd(dev_niak, incr);
         pos = __shfl_sync(ALL_LANES, pos, 0);
         if (pos + incr <= nak * Spatial2::LSTCAP) {
            if (ilane < incr)
               iak[pos + ilane] = wy;
            for (int i = 0; i < incr; ++i) {
               int bufp = i * WARP_SIZE + ilane;
               int num = bufp < nknb ? buffer[bufp] : 0;
               lst[(pos + i) * WARP_SIZE + ilane] = num;
            }
         }
      }
   }
}


struct spatial2_less
{
   __device__
   bool operator()(int2 a, int2 b)
   {
      return a.x < b.x;
   }
};


template <class IMG>
void run_spatial2_step5(Spatial2Unit u)
{
   u->niak = 0;
   int* dev_niak = &u->update[2];
   real cutbuf = u->cutoff + u->buffer;
   launch_k1s(g::s0, u->nakp, spatial2_step5<IMG>, //
              u->bnum, u->iakpl_rev, u->nstype, u->si1, u->si2, u->si3,
              u->si4,                               //
              dev_niak, u->iak, u->lst,             //
              n, u->nak, cutbuf, TINKER_IMAGE_ARGS, //
              u->akpf, u->sorted, u->akc, u->half);
   darray::copyout(g::q0, 1, &u->niak, dev_niak);
   wait_for(g::q0);
   if (u->niak > u->nak * Spatial2::LSTCAP) {
      int cap = Spatial2::LSTCAP;
      TINKER_THROW(
         format("An internal array in Spatial2 requested %1$d elements, "
                "but only %4$d (%3$d*%2$d) were allocated. "
                "Please increase Spatial2::LSTCAP (current value %3$d) "
                "so as to make Spatial2::LSTCAP*%2$d >= %1$d.\n",
                u->niak, u->nak, cap, cap * u->nak));
   }
}


void spatial_data_init_cu(Spatial2Unit u)
{
   u->fresh = -1; // 0xFFFFFFFF;
   const real lbuf = u->buffer;
   const int n = u->n;


   auto policy = thrust::cuda::par(thrust_cache).on(g::s0);


   const auto* lx = u->x;
   const auto* ly = u->y;
   const auto* lz = u->z;
   int ZERO_LBUF = (lbuf <= 0 ? 1 : 0);
   real cutbuf = (u->cutoff + lbuf);
   int2* b2num = (int2*)u->update;
   launch_k1s(g::s0, n, spatial2_step1, //
              n, u->pz, b2num, lx, ly, lz, TINKER_IMAGE_ARGS, u->nakpk,
              u->akpf);
   // thrust::sort(policy, b2num, b2num + n, spatial2_less());
   thrust::stable_sort(policy, b2num, b2num + n, spatial2_less());
   launch_k1s(g::s0, n, spatial2_step2, //
              n, u->sorted, u->bnum, b2num, lx, ly, lz, ZERO_LBUF, u->xold,
              u->yold, u->zold, //
              TINKER_IMAGE_ARGS, cutbuf, u->akc, u->half);


   auto& si1 = u->si1;
   auto& si2 = u->si2;
   auto& si3 = u->si3;
   auto& si4 = u->si4;
   int* nakpl_ptr0 = &u->update[0];
   launch_k1s(g::s0, n, spatial2_step3,       //
              u->nak, u->akpf, nakpl_ptr0,    //
              u->bnum, u->nstype,             //
              si1.ns, si1.js, si2.ns, si2.js, //
              si3.ns, si3.js, si4.ns, si4.js);
   darray::copyout(g::q0, 1, &u->nakpl, nakpl_ptr0);
   wait_for(g::q0);
   if (WARP_SIZE + u->nakpl > u->cap_nakpl) {
      u->cap_nakpl = WARP_SIZE + u->nakpl;
      darray::deallocate(u->iakpl);
      darray::allocate(u->cap_nakpl, &u->iakpl);
      if (u->nstype >= 1) {
         auto& si = si1;
         darray::deallocate(si.bit0);
         darray::allocate(32 * u->cap_nakpl, &si.bit0);
      }
      if (u->nstype >= 2) {
         auto& si = si2;
         darray::deallocate(si.bit0);
         darray::allocate(32 * u->cap_nakpl, &si.bit0);
      }
      if (u->nstype >= 3) {
         auto& si = si3;
         darray::deallocate(si.bit0);
         darray::allocate(32 * u->cap_nakpl, &si.bit0);
      }
      if (u->nstype >= 4) {
         auto& si = si4;
         darray::deallocate(si.bit0);
         darray::allocate(32 * u->cap_nakpl, &si.bit0);
      }
   }


   int* nakpl_ptr1 = &u->update[1];
   launch_k1s(g::s0, u->nakp, spatial2_step4,                        //
              u->nakpk, nakpl_ptr1, u->akpf, u->iakpl, u->iakpl_rev, //
              u->cap_nakpl, u->nstype,                               //
              si1.bit0, si2.bit0, si3.bit0, si4.bit0);


   if (box_shape == ORTHO_BOX) {
      run_spatial2_step5<PBC_ORTHO>(u);
   } else if (box_shape == MONO_BOX) {
      run_spatial2_step5<PBC_MONO>(u);
   } else if (box_shape == TRI_BOX) {
      run_spatial2_step5<PBC_TRI>(u);
   } else if (box_shape == OCT_BOX) {
      run_spatial2_step5<PBC_OCT>(u);
   } else {
      assert(false);
   }
}


void spatial_data_update_sorted(Spatial2Unit u)
{
   launch_k1s(g::s0, u->n, spatial2_update_sorted, //
              u->n, u->sorted, u->x, u->y, u->z,   //
              TINKER_IMAGE_ARGS, u->cutoff, u->akc, u->half);
}
}
