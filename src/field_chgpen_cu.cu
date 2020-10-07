#include "add.h"
#include "empole_chgpen.h"
#include "epolar_chgpen.h"
#include "glob.spatial.h"
#include "image.h"
#include "launch.h"
#include "md.h"
#include "pme.h"
#include "seq_pair_field_chgpen.h"
#include "seq_triangle.h"
#include "switch.h"
#include "tool/cudalib.h"
#include "tool/gpu_card.h"


namespace tinker {
// #define DFIELDPARAS                                                            \
//    real(*restrict field)[3], real *restrict pcore, real *restrict pval,        \
//       real *restrict palpha, const real(*restrict rpole)[10],                  \
//       TINKER_IMAGE_PARAMS, real off2


// template <class ETYP>
// __launch_bounds__(BLOCK_DIM) __global__
// void dfield_chgpen_cu1(DFIELDPARAS, int n,
//                        const Spatial::SortedAtom* restrict sorted, int niak,
//                        const int* restrict iak, const int* restrict lst,
//                        real aewald)
// {
//    const int iwarp = (threadIdx.x + blockIdx.x * blockDim.x) / WARP_SIZE;
//    const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
//    const int ilane = threadIdx.x & (WARP_SIZE - 1);


//    struct Data
//    {
//       real3 fkd;
//       real3 rk;
//       real ck;
//       real dkx, dky, dkz;
//       real qkxx, qkxy, qkxz, qkyy, qkyz, qkzz;
//       real core, val, alpha;
//    };
//    __shared__ Data data[BLOCK_DIM];


//    for (int iw = iwarp; iw < niak; iw += nwarp) {
//       real3 fid = make_real3(0, 0, 0);
//       int atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
//       real3 ri = make_real3(sorted[atomi].x, sorted[atomi].y,
//       sorted[atomi].z); int i = sorted[atomi].unsorted; real ci =
//       rpole[i][mpl_pme_0]; real dix = rpole[i][mpl_pme_x]; real diy =
//       rpole[i][mpl_pme_y]; real diz = rpole[i][mpl_pme_z]; real qixx =
//       rpole[i][mpl_pme_xx]; real qixy = rpole[i][mpl_pme_xy]; real qixz =
//       rpole[i][mpl_pme_xz]; real qiyy = rpole[i][mpl_pme_yy]; real qiyz =
//       rpole[i][mpl_pme_yz]; real qizz = rpole[i][mpl_pme_zz]; real corei =
//       pcore[i]; real alphai = palpha[i]; real vali = pval[i];


//       data[threadIdx.x].fkd = make_real3(0, 0, 0);
//       int shatomk = lst[iw * WARP_SIZE + ilane];
//       data[threadIdx.x].rk =
//          make_real3(sorted[shatomk].x, sorted[shatomk].y, sorted[shatomk].z);
//       int shk = sorted[shatomk].unsorted;
//       data[threadIdx.x].ck = rpole[shk][mpl_pme_0];
//       data[threadIdx.x].dkx = rpole[shk][mpl_pme_x];
//       data[threadIdx.x].dky = rpole[shk][mpl_pme_y];
//       data[threadIdx.x].dkz = rpole[shk][mpl_pme_z];
//       data[threadIdx.x].qkxx = rpole[shk][mpl_pme_xx];
//       data[threadIdx.x].qkxy = rpole[shk][mpl_pme_xy];
//       data[threadIdx.x].qkxz = rpole[shk][mpl_pme_xz];
//       data[threadIdx.x].qkyy = rpole[shk][mpl_pme_yy];
//       data[threadIdx.x].qkyz = rpole[shk][mpl_pme_yz];
//       data[threadIdx.x].qkzz = rpole[shk][mpl_pme_zz];
//       data[threadIdx.x].core = pcore[shk];
//       data[threadIdx.x].val = pval[shk];
//       data[threadIdx.x].alpha = palpha[shk];


//       for (int j = 0; j < WARP_SIZE; ++j) {
//          int srclane = (ilane + j) & (WARP_SIZE - 1);
//          int klane = srclane + threadIdx.x - ilane;
//          int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
//          real3 dr = data[klane].rk - ri;


//          real r2 = image2(dr.x, dr.y, dr.z);
//          if (atomi < atomk && r2 <= off2) {
//             if CONSTEXPR (eq<ETYP, EWALD>()) {
//                pair_dfield_chgpen<EWALD>(
//                   r2, dr.x, dr.y, dr.z, 1, ci, dix, diy, diz, corei, vali,
//                   alphai, qixx, qixy, qixz, qiyy, qiyz, qizz, data[klane].ck,
//                   data[klane].dkx, data[klane].dky, data[klane].dkz,
//                   data[klane].core, data[klane].val, data[klane].alpha,
//                   data[klane].qkxx, data[klane].qkxy, data[klane].qkxz,
//                   data[klane].qkyy, data[klane].qkyz, data[klane].qkzz,
//                   aewald, fid, data[klane].fkd);
//             }
//             if CONSTEXPR (eq<ETYP, NON_EWALD>()) {
//                pair_dfield_chgpen<NON_EWALD>(
//                   r2, dr.x, dr.y, dr.z, 1, ci, dix, diy, diz, corei, vali,
//                   alphai, qixx, qixy, qixz, qiyy, qiyz, qizz, data[klane].ck,
//                   data[klane].dkx, data[klane].dky, data[klane].dkz,
//                   data[klane].core, data[klane].val, data[klane].alpha,
//                   data[klane].qkxx, data[klane].qkxy, data[klane].qkxz,
//                   data[klane].qkyy, data[klane].qkyz, data[klane].qkzz, 0,
//                   fid, data[klane].fkd);
//             }
//          } // end if (include)
//       }


//       atomic_add(fid.x, &field[i][0]);
//       atomic_add(fid.y, &field[i][1]);
//       atomic_add(fid.z, &field[i][2]);
//       atomic_add(data[threadIdx.x].fkd.x, &field[shk][0]);
//       atomic_add(data[threadIdx.x].fkd.y, &field[shk][1]);
//       atomic_add(data[threadIdx.x].fkd.z, &field[shk][2]);
//    } // end for (iw)
// }


// __global__
// void dfield_chgpen_cu2(DFIELDPARAS, const real* restrict x,
//                        const real* restrict y, const real* restrict z,
//                        int ndexclude, const int (*restrict dexclude)[2],
//                        const real(*restrict dexclude_scale))
// {
//    for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < ndexclude;
//         ii += blockDim.x * gridDim.x) {
//       int i = dexclude[ii][0];
//       int k = dexclude[ii][1];
//       real dscale = dexclude_scale[ii];

//       real xi = x[i];
//       real yi = y[i];
//       real zi = z[i];
//       real ci = rpole[i][mpl_pme_0];
//       real dix = rpole[i][mpl_pme_x];
//       real diy = rpole[i][mpl_pme_y];
//       real diz = rpole[i][mpl_pme_z];
//       real qixx = rpole[i][mpl_pme_xx];
//       real qixy = rpole[i][mpl_pme_xy];
//       real qixz = rpole[i][mpl_pme_xz];
//       real qiyy = rpole[i][mpl_pme_yy];
//       real qiyz = rpole[i][mpl_pme_yz];
//       real qizz = rpole[i][mpl_pme_zz];
//       real corei = pcore[i];
//       real alphai = palpha[i];
//       real vali = pval[i];


//       real xr = x[k] - xi;
//       real yr = y[k] - yi;
//       real zr = z[k] - zi;


//       real r2 = image2(xr, yr, zr);
//       if (r2 <= off2) {
//          real3 fid = make_real3(0, 0, 0);
//          real3 fkd = make_real3(0, 0, 0);
//          pair_dfield_chgpen<NON_EWALD>(
//             r2, xr, yr, zr, dscale, ci, dix, diy, diz, corei, vali, alphai,
//             qixx, qixy, qixz, qiyy, qiyz, qizz, rpole[k][mpl_pme_0],
//             rpole[k][mpl_pme_x], rpole[k][mpl_pme_y], rpole[k][mpl_pme_z],
//             pcore[k], pval[k], palpha[k], rpole[k][mpl_pme_xx],
//             rpole[k][mpl_pme_xy], rpole[k][mpl_pme_xz], rpole[k][mpl_pme_yy],
//             rpole[k][mpl_pme_yz], rpole[k][mpl_pme_zz], 0, fid, fkd);


//          atomic_add(fid.x, &field[i][0]);
//          atomic_add(fid.y, &field[i][1]);
//          atomic_add(fid.z, &field[i][2]);


//          atomic_add(fkd.x, &field[k][0]);
//          atomic_add(fkd.y, &field[k][1]);
//          atomic_add(fkd.z, &field[k][2]);
//       } // end if (include)
//    }
// }

template <class ETYP>
__launch_bounds__(BLOCK_DIM) __global__
void dfield_chgpen_cu1(
   int n, TINKER_IMAGE_PARAMS


   ,
   real off, const unsigned* restrict dinfo, int nexclude,
   const int (*restrict exclude)[2], const real (*restrict exclude_scale)[3],
   const real* restrict x, const real* restrict y, const real* restrict z,
   const Spatial::SortedAtom* restrict sorted, int nakpl,
   const int* restrict iakpl, int niak, const int* restrict iak,
   const int* restrict lst, real (*restrict field)[3],
   const real (*restrict rpole)[10], real* restrict pcore, real* restrict pval,
   const real* restrict palpha, real aewald)
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
   __shared__ real shcorei[BLOCK_DIM];
   __shared__ real shalphai[BLOCK_DIM];
   __shared__ real shvali[BLOCK_DIM];
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
   real corek;
   real alphak;
   real valk;


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
      real scaleb = exclude_scale[ii][1];


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
      real corei = pcore[shi];
      real alphai = palpha[shi];
      real vali = pval[shi];
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
      corek = pcore[k];
      alphak = palpha[k];
      valk = pval[k];


      int klane = threadIdx.x;
      constexpr bool incl = true;
      real xr = xk - xi;
      real yr = yk - yi;
      real zr = zk - zi;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         pair_dfield_chgpen<ETYP>(
            r2, xr, yr, zr, scaleb, ci, dix, diy, diz, corei, vali, alphai,
            qixx, qixy, qixz, qiyy, qiyz, qizz, ck, dkx, dky, dkz, corek, valk,
            alphak, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, aewald, shfidx[klane],
            shfidy[klane], shfidz[klane], fkdx, fkdy, fkdz);
      } // end if (include)


      atomic_add(shfidx[threadIdx.x], &field[shi][0]);
      atomic_add(shfidy[threadIdx.x], &field[shi][1]);
      atomic_add(shfidz[threadIdx.x], &field[shi][2]);
      atomic_add(fkdx, &field[k][0]);
      atomic_add(fkdy, &field[k][1]);
      atomic_add(fkdz, &field[k][2]);
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
      shcorei[threadIdx.x] = pcore[shi];
      shalphai[threadIdx.x] = palpha[shi];
      shvali[threadIdx.x] = pval[shi];
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
      corek = pcore[k];
      alphak = palpha[k];
      valk = pval[k];


      unsigned int dinfo0 = dinfo[iw * WARP_SIZE + ilane];


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
         real corei = shcorei[klane];
         real alphai = shalphai[klane];
         real vali = shvali[klane];


         bool incl = iid < kid and kid < n;
         incl = incl and (dinfo0 & srcmask) == 0;
         real scaleb = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_dfield_chgpen<ETYP>(
               r2, xr, yr, zr, scaleb, ci, dix, diy, diz, corei, vali, alphai,
               qixx, qixy, qixz, qiyy, qiyz, qizz, ck, dkx, dky, dkz, corek,
               valk, alphak, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, aewald,
               shfidx[klane], shfidy[klane], shfidz[klane], fkdx, fkdy, fkdz);
         } // end if (include)


         shiid = __shfl_sync(ALL_LANES, shiid, ilane + 1);
      }


      atomic_add(shfidx[threadIdx.x], &field[shi][0]);
      atomic_add(shfidy[threadIdx.x], &field[shi][1]);
      atomic_add(shfidz[threadIdx.x], &field[shi][2]);
      atomic_add(fkdx, &field[k][0]);
      atomic_add(fkdy, &field[k][1]);
      atomic_add(fkdz, &field[k][2]);
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
      shcorei[threadIdx.x] = pcore[shi];
      shalphai[threadIdx.x] = palpha[shi];
      shvali[threadIdx.x] = pval[shi];
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
      corek = pcore[k];
      alphak = palpha[k];
      valk = pval[k];


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
         real corei = shcorei[klane];
         real alphai = shalphai[klane];
         real vali = shvali[klane];


         bool incl = atomk > 0;
         real scaleb = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_dfield_chgpen<ETYP>(
               r2, xr, yr, zr, scaleb, ci, dix, diy, diz, corei, vali, alphai,
               qixx, qixy, qixz, qiyy, qiyz, qizz, ck, dkx, dky, dkz, corek,
               valk, alphak, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, aewald,
               shfidx[klane], shfidy[klane], shfidz[klane], fkdx, fkdy, fkdz);
         } // end if (include)
      }


      atomic_add(shfidx[threadIdx.x], &field[shi][0]);
      atomic_add(shfidy[threadIdx.x], &field[shi][1]);
      atomic_add(shfidz[threadIdx.x], &field[shi][2]);
      atomic_add(fkdx, &field[k][0]);
      atomic_add(fkdy, &field[k][1]);
      atomic_add(fkdz, &field[k][2]);
   }
   // */
} // generated by ComplexKernelBuilder (ck.py) 1.5.1

template <class ETYP>
void dfield_chgpen_cu(real (*field)[3])
{
   // const auto& st = *mspatial_unit;
   // const real off = switch_off(switch_ewald);
   // const real off2 = off * off;  

   const auto& st = *mspatial_v2_unit;
   real off;

   if CONSTEXPR (eq<ETYP, EWALD>())
      off = switch_off(switch_ewald);
   else
      off = switch_off(switch_mpole);


   real aewald = 0;
   if CONSTEXPR (eq<ETYP, EWALD>()) {
      PMEUnit pu = ppme_unit;
      aewald = pu->aewald;
   }

   // if (st.niak > 0) {
   //    launch_k1s(nonblk, WARP_SIZE * st.niak, dfield_chgpen_cu1<EWALD>,
   //    field,
   //               pcore, pval, palpha, rpole, TINKER_IMAGE_ARGS, off2, n,
   //               st.sorted, st.niak, st.iak, st.lst, aewald);
   // }
   // if (ndexclude > 0) {
   //    launch_k1s(nonblk, ndexclude, dfield_chgpen_cu2, field, pcore, pval,
   //               palpha, rpole, TINKER_IMAGE_ARGS, off2, x, y, z, ndexclude,
   //               dexclude, dexclude_scale);
   // }
   int ngrid = get_grid_size(BLOCK_DIM);
   dfield_chgpen_cu1<ETYP><<<ngrid, BLOCK_DIM, 0, nonblk>>>(
      st.n, TINKER_IMAGE_ARGS, off, st.si1.bit0, nmdwexclude, mdwexclude,
      mdwexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl,
      st.niak, st.iak, st.lst, field, rpole, pcore, pval, palpha, aewald);
}

void dfield_chgpen_ewald_real_cu(real (*field)[3]) {
   dfield_chgpen_cu<EWALD>(field);
}
void dfield_chgpen_nonewald_cu(real (*field)[3])
{
   darray::zero(PROCEED_NEW_Q, n, field);

   dfield_chgpen_cu<NON_EWALD>(field);
   // const auto& st = *mspatial_unit;
   // const real off = switch_off(switch_mpole);
   // const real off2 = off * off;
   // const auto& st = *mspatial_v2_unit;
   // real off;
   
   // if (st.niak > 0) {
   //    launch_k1s(nonblk, WARP_SIZE * st.niak, dfield_chgpen_cu1<NON_EWALD>,
   //               field, pcore, pval, palpha, rpole, TINKER_IMAGE_ARGS, off2,
   //               n, st.sorted, st.niak, st.iak, st.lst, 0);
   // }
   // if (ndexclude > 0) {
   //    launch_k1s(nonblk, ndexclude, dfield_chgpen_cu2, field, pcore, pval,
   //               palpha, rpole, TINKER_IMAGE_ARGS, off2, x, y, z, ndexclude,
   //               dexclude, dexclude_scale);
   // }

   // int ngrid = get_grid_size(BLOCK_DIM);
   // dfield_chgpen_cu1<NON_EWALD><<<ngrid, BLOCK_DIM, 0, nonblk>>>(
   //    st.n, TINKER_IMAGE_ARGS, off, st.si1.bit0, nmdwexclude, mdwexclude,
   //    mdwexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl,
   //    st.niak, st.iak, st.lst, field, rpole, pcore, pval, palpha, 0);
}


// #define UFIELDPARAS                                                            \
//    const real(*restrict uind)[3], real(*restrict field)[3],                    \
//       real *restrict pcore, real *restrict pval, real *restrict palpha,        \
//       TINKER_IMAGE_PARAMS, real off2


// template <class ETYP>
// __launch_bounds__(BLOCK_DIM) __global__
// void ufield_chgpen_cu1(UFIELDPARAS, int n,
//                        const Spatial::SortedAtom* restrict sorted, int niak,
//                        const int* restrict iak, const int* restrict lst,
//                        real aewald)
// {
//    const int iwarp = (threadIdx.x + blockIdx.x * blockDim.x) / WARP_SIZE;
//    const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
//    const int ilane = threadIdx.x & (WARP_SIZE - 1);


//    struct Data
//    {
//       real3 fkd, ukd, rk;
//       real core, alpha, val;
//    };
//    __shared__ Data data[BLOCK_DIM];


//    for (int iw = iwarp; iw < niak; iw += nwarp) {
//       real3 fid = make_real3(0, 0, 0);
//       int atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
//       real3 ri = make_real3(sorted[atomi].x, sorted[atomi].y,
//       sorted[atomi].z); int i = sorted[atomi].unsorted; real3 uid =
//       make_real3(uind[i][0], uind[i][1], uind[i][2]); real corei = pcore[i];
//       real alphai = palpha[i];
//       real vali = pval[i];


//       data[threadIdx.x].fkd = make_real3(0, 0, 0);
//       int shatomk = lst[iw * WARP_SIZE + ilane];
//       data[threadIdx.x].rk =
//          make_real3(sorted[shatomk].x, sorted[shatomk].y, sorted[shatomk].z);
//       int shk = sorted[shatomk].unsorted;
//       data[threadIdx.x].ukd =
//          make_real3(uind[shk][0], uind[shk][1], uind[shk][2]);
//       data[threadIdx.x].core = pcore[shk];
//       data[threadIdx.x].val = pval[shk];
//       data[threadIdx.x].alpha = palpha[shk];


//       for (int j = 0; j < WARP_SIZE; ++j) {
//          int srclane = (ilane + j) & (WARP_SIZE - 1);
//          int klane = srclane + threadIdx.x - ilane;
//          int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
//          real3 dr = data[klane].rk - ri;


//          real r2 = image2(dr.x, dr.y, dr.z);
//          if (atomi < atomk && r2 <= off2) {
//             if CONSTEXPR (eq<ETYP, EWALD>()) {
//                pair_ufield_chgpen<EWALD>(
//                   r2, dr.x, dr.y, dr.z, 1, uid.x, uid.y, uid.z, corei, vali,
//                   alphai, data[klane].ukd.x, data[klane].ukd.y,
//                   data[klane].ukd.z, data[klane].core, data[klane].val,
//                   data[klane].alpha, aewald, fid, data[klane].fkd);
//             }
//             if CONSTEXPR (eq<ETYP, NON_EWALD>()) {
//                pair_ufield_chgpen<NON_EWALD>(
//                   r2, dr.x, dr.y, dr.z, 1, uid.x, uid.y, uid.z, corei, vali,
//                   alphai, data[klane].ukd.x, data[klane].ukd.y,
//                   data[klane].ukd.z, data[klane].core, data[klane].val,
//                   data[klane].alpha, aewald, fid, data[klane].fkd);
//             }
//          } // end if (include)
//       }


//       atomic_add(fid.x, &field[i][0]);
//       atomic_add(fid.y, &field[i][1]);
//       atomic_add(fid.z, &field[i][2]);
//       atomic_add(data[threadIdx.x].fkd.x, &field[shk][0]);
//       atomic_add(data[threadIdx.x].fkd.y, &field[shk][1]);
//       atomic_add(data[threadIdx.x].fkd.z, &field[shk][2]);
//    } // end for (iw)
// }


// __global__
// void ufield_chgpen_cu2(UFIELDPARAS, const real* restrict x,
//                        const real* restrict y, const real* restrict z,
//                        int ndwexclude, const int (*restrict dwexclude)[2],
//                        const real* restrict dwexclude_scale)
// {
//    for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < ndwexclude;
//         ii += blockDim.x * gridDim.x) {
//       int i = dwexclude[ii][0];
//       int k = dwexclude[ii][1];
//       real wscale = dwexclude_scale[ii] - 1;

//       real xi = x[i];
//       real yi = y[i];
//       real zi = z[i];
//       real3 uid = make_real3(uind[i][0], uind[i][1], uind[i][2]);
//       real corei = pcore[i];
//       real alphai = palpha[i];
//       real vali = pval[i];


//       real xr = x[k] - xi;
//       real yr = y[k] - yi;
//       real zr = z[k] - zi;


//       real r2 = image2(xr, yr, zr);
//       if (r2 <= off2) {
//          real3 fid = make_real3(0, 0, 0);
//          real3 fkd = make_real3(0, 0, 0);
//          pair_ufield_chgpen<NON_EWALD>(r2, xr, yr, zr, wscale, uid.x, uid.y,
//                                        uid.z, corei, vali, alphai,
//                                        uind[k][0], uind[k][1], uind[k][2],
//                                        pcore[k], pval[k], palpha[k], 0, fid,
//                                        fkd);


//          atomic_add(fid.x, &field[i][0]);
//          atomic_add(fid.y, &field[i][1]);
//          atomic_add(fid.z, &field[i][2]);


//          atomic_add(fkd.x, &field[k][0]);
//          atomic_add(fkd.y, &field[k][1]);
//          atomic_add(fkd.z, &field[k][2]);
//       } // end if (include)
//    }
// }
template <class ETYP>
__launch_bounds__(BLOCK_DIM) __global__
void ufield_chgpen_cu1(
   int n, TINKER_IMAGE_PARAMS,
   real off, const unsigned* restrict winfo, int nexclude,
   const int (*restrict exclude)[2], const real (*restrict exclude_scale)[3],
   const real* restrict x, const real* restrict y, const real* restrict z,
   const Spatial::SortedAtom* restrict sorted, int nakpl,
   const int* restrict iakpl, int niak, const int* restrict iak,
   const int* restrict lst, const real (*restrict uind)[3],
   real (*restrict field)[3], real* restrict pcore, real* restrict pval,
   const real* restrict palpha, real aewald)
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
   __shared__ real shcorei[BLOCK_DIM];
   __shared__ real shalphai[BLOCK_DIM];
   __shared__ real shvali[BLOCK_DIM];
   real ukdx;
   real ukdy;
   real ukdz;
   real corek;
   real alphak;
   real valk;


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
      real scalec = exclude_scale[ii][2];


      real xi = x[shi];
      real yi = y[shi];
      real zi = z[shi];
      xk = x[k];
      yk = y[k];
      zk = z[k];
      real uidx = uind[shi][0];
      real uidy = uind[shi][1];
      real uidz = uind[shi][2];
      real corei = pcore[shi];
      real alphai = palpha[shi];
      real vali = pval[shi];
      ukdx = uind[k][0];
      ukdy = uind[k][1];
      ukdz = uind[k][2];
      corek = pcore[k];
      alphak = palpha[k];
      valk = pval[k];


      int klane = threadIdx.x;
      constexpr bool incl = true;
      real xr = xk - xi;
      real yr = yk - yi;
      real zr = zk - zi;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         pair_ufield_chgpen<ETYP>(
            r2, xr, yr, zr, scalec, uidx, uidy, uidz, corei, vali, alphai, ukdx,
            ukdy, ukdz, corek, valk, alphak, aewald, shfidx[klane],
            shfidy[klane], shfidz[klane], fkdx, fkdy, fkdz);
      } // end if (include)


      atomic_add(shfidx[threadIdx.x], &field[shi][0]);
      atomic_add(shfidy[threadIdx.x], &field[shi][1]);
      atomic_add(shfidz[threadIdx.x], &field[shi][2]);
      atomic_add(fkdx, &field[k][0]);
      atomic_add(fkdy, &field[k][1]);
      atomic_add(fkdz, &field[k][2]);
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


      shuidx[threadIdx.x] = uind[shi][0];
      shuidy[threadIdx.x] = uind[shi][1];
      shuidz[threadIdx.x] = uind[shi][2];
      shcorei[threadIdx.x] = pcore[shi];
      shalphai[threadIdx.x] = palpha[shi];
      shvali[threadIdx.x] = pval[shi];
      ukdx = uind[k][0];
      ukdy = uind[k][1];
      ukdz = uind[k][2];
      corek = pcore[k];
      alphak = palpha[k];
      valk = pval[k];


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
         real corei = shcorei[klane];
         real alphai = shalphai[klane];
         real vali = shvali[klane];


         bool incl = iid < kid and kid < n;
         incl = incl and (winfo0 & srcmask) == 0;
         real scalec = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_ufield_chgpen<ETYP>(
               r2, xr, yr, zr, scalec, uidx, uidy, uidz, corei, vali, alphai,
               ukdx, ukdy, ukdz, corek, valk, alphak, aewald, shfidx[klane],
               shfidy[klane], shfidz[klane], fkdx, fkdy, fkdz);
         } // end if (include)


         shiid = __shfl_sync(ALL_LANES, shiid, ilane + 1);
      }


      atomic_add(shfidx[threadIdx.x], &field[shi][0]);
      atomic_add(shfidy[threadIdx.x], &field[shi][1]);
      atomic_add(shfidz[threadIdx.x], &field[shi][2]);
      atomic_add(fkdx, &field[k][0]);
      atomic_add(fkdy, &field[k][1]);
      atomic_add(fkdz, &field[k][2]);
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


      shuidx[threadIdx.x] = uind[shi][0];
      shuidy[threadIdx.x] = uind[shi][1];
      shuidz[threadIdx.x] = uind[shi][2];
      shcorei[threadIdx.x] = pcore[shi];
      shalphai[threadIdx.x] = palpha[shi];
      shvali[threadIdx.x] = pval[shi];
      ukdx = uind[k][0];
      ukdy = uind[k][1];
      ukdz = uind[k][2];
      corek = pcore[k];
      alphak = palpha[k];
      valk = pval[k];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         real xi = shxi[klane];
         real yi = shyi[klane];
         real zi = shzi[klane];
         real uidx = shuidx[klane];
         real uidy = shuidy[klane];
         real uidz = shuidz[klane];
         real corei = shcorei[klane];
         real alphai = shalphai[klane];
         real vali = shvali[klane];


         bool incl = atomk > 0;
         real scalec = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;
         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_ufield_chgpen<ETYP>(
               r2, xr, yr, zr, scalec, uidx, uidy, uidz, corei, vali, alphai,
               ukdx, ukdy, ukdz, corek, valk, alphak, aewald, shfidx[klane],
               shfidy[klane], shfidz[klane], fkdx, fkdy, fkdz);
         } // end if (include)
      }


      atomic_add(shfidx[threadIdx.x], &field[shi][0]);
      atomic_add(shfidy[threadIdx.x], &field[shi][1]);
      atomic_add(shfidz[threadIdx.x], &field[shi][2]);
      atomic_add(fkdx, &field[k][0]);
      atomic_add(fkdy, &field[k][1]);
      atomic_add(fkdz, &field[k][2]);
   }
   // */
} // generated by ComplexKernelBuilder (ck.py) 1.5.1


template <class ETYP>
void ufield_chgpen_cu(const real (*uind)[3], real (*field)[3]) {
   
   const auto& st = *mspatial_v2_unit;
   real off;

   if CONSTEXPR (eq<ETYP, EWALD>())
      off = switch_off(switch_ewald);
   else
      off = switch_off(switch_mpole);


   real aewald = 0;
   if CONSTEXPR (eq<ETYP, EWALD>()) {
      PMEUnit pu = ppme_unit;
      aewald = pu->aewald;
   }

   int ngrid = get_grid_size(BLOCK_DIM);
   ufield_chgpen_cu1<ETYP><<<ngrid, BLOCK_DIM, 0, nonblk>>>(
      st.n, TINKER_IMAGE_ARGS, off, st.si1.bit0, nmdwexclude, mdwexclude,
      mdwexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl,
      st.niak, st.iak, st.lst, uind, field, pcore, pval, palpha, aewald);
}

void ufield_chgpen_ewald_real_cu(const real (*uind)[3], real (*field)[3])
{
   ufield_chgpen_cu<EWALD>(uind,field);
   // const auto& st = *mspatial_unit;
   // const real off = switch_off(switch_ewald);
   // const real off2 = off * off;

   // const PMEUnit pu = ppme_unit;
   // const real aewald = pu->aewald;


   // if (st.niak > 0) {
   //    launch_k1s(nonblk, WARP_SIZE * st.niak, ufield_chgpen_cu1<EWALD>, uind,
   //               field, pcore, pval, palpha, TINKER_IMAGE_ARGS, off2, n,
   //               st.sorted, st.niak, st.iak, st.lst, aewald);
   // }
   // if (ndwexclude) {
   //    launch_k1s(nonblk, ndwexclude, ufield_chgpen_cu2, uind, field, pcore,
   //    pval,
   //               palpha, TINKER_IMAGE_ARGS, off2, x, y, z, ndwexclude,
   //               dwexclude, dwexclude_scale);
   
}


void ufield_chgpen_nonewald_cu(const real (*uind)[3], real (*field)[3])
{
   darray::zero(PROCEED_NEW_Q, n, field);
   ufield_chgpen_cu<NON_EWALD>(uind,field);

   // const auto& st = *mspatial_unit;
   // const real off = switch_off(switch_mpole);
   // const real off2 = off * off;


   // if (st.niak > 0) {
   //    launch_k1s(nonblk, WARP_SIZE * st.niak, ufield_chgpen_cu1<NON_EWALD>,
   //               uind, field, pcore, pval, palpha, TINKER_IMAGE_ARGS, off2,
   //               n, st.sorted, st.niak, st.iak, st.lst, 0);
   // }
   // if (ndwexclude) {
   //    launch_k1s(nonblk, ndwexclude, ufield_chgpen_cu2, uind, field, pcore,
   //    pval,
   //               palpha, TINKER_IMAGE_ARGS, off2, x, y, z, ndwexclude,
   //               dwexclude, dwexclude_scale);

}
}
