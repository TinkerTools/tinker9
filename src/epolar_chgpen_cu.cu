#include "add.h"
#include "empole_chgpen.h"
#include "epolar_chgpen.h"
#include "epolar_trq.h"
#include "glob.spatial.h"
#include "image.h"
#include "launch.h"
#include "md.h"
#include "pmestuf.h"
#include "seq_bsplgen.h"
#include "seq_pair_polar_chgpen.h"
#include "seq_triangle.h"
#include "switch.h"
#include "tool/gpu_card.h"


namespace tinker {
// #define POLARPARAS                                                             \
//    size_t bufsize, count_buffer restrict nep, energy_buffer restrict ep,       \
//       virial_buffer restrict vir_ep, grad_prec *restrict gx,                   \
//       grad_prec *restrict gy, grad_prec *restrict gz, real(*restrict ufld)[3], \
//       real(*restrict dufld)[6], TINKER_IMAGE_PARAMS, real off2, real f,        \
//       const real(*restrict rpole)[10], real *restrict pcore,                   \
//       real *restrict pval, real *restrict palpha,                              \
//       const real(*restrict uind)[3]


// template <class Ver, class ETYP>
// __global__
// void epolar_chgpen_cu1(POLARPARAS, const Spatial::SortedAtom* restrict
// sorted,
//                        int niak, const int* restrict iak,
//                        const int* restrict lst, int n, real aewald)
// {
//    constexpr bool do_e = Ver::e;
//    constexpr bool do_a = Ver::a;
//    constexpr bool do_g = Ver::g;
//    constexpr bool do_v = Ver::v;


//    const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
//    const int iwarp = ithread / WARP_SIZE;
//    const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
//    const int ilane = threadIdx.x & (WARP_SIZE - 1);
//    const int offset = ithread & (bufsize - 1);


//    MAYBE_UNUSED int ctl;
//    MAYBE_UNUSED real etl;
//    MAYBE_UNUSED real vtlxx, vtlxy, vtlxz, vtlyy, vtlyz, vtlzz;
//    MAYBE_UNUSED real gxi, gyi, gzi, txi, tyi, tzi, dui[6];
//    MAYBE_UNUSED __shared__ real gxk[BLOCK_DIM], gyk[BLOCK_DIM],
//    gzk[BLOCK_DIM],
//       txk[BLOCK_DIM], tyk[BLOCK_DIM], tzk[BLOCK_DIM], duk[BLOCK_DIM][6];


//    for (int iw = iwarp; iw < niak; iw += nwarp) {
//       if CONSTEXPR (do_a)
//          ctl = 0;
//       if CONSTEXPR (do_e)
//          etl = 0;

//       if CONSTEXPR (do_v) {
//          vtlxx = 0;
//          vtlxy = 0;
//          vtlxz = 0;
//          vtlyy = 0;
//          vtlyz = 0;
//          vtlzz = 0;
//       }
//       if CONSTEXPR (do_g) {
//          gxi = 0;
//          gyi = 0;
//          gzi = 0;
//          txi = 0;
//          tyi = 0;
//          tzi = 0;
//          #pragma unroll
//          for (int i = 0; i < 6; ++i) {
//             dui[i] = 0;
//          }
//          gxk[threadIdx.x] = 0;
//          gyk[threadIdx.x] = 0;
//          gzk[threadIdx.x] = 0;
//          txk[threadIdx.x] = 0;
//          tyk[threadIdx.x] = 0;
//          tzk[threadIdx.x] = 0;
//          #pragma unroll
//          for (int i = 0; i < 6; ++i) {
//             duk[threadIdx.x][i] = 0;
//          }
//       }


//       int atomi;
//       atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
//       real xi = sorted[atomi].x;
//       real yi = sorted[atomi].y;
//       real zi = sorted[atomi].z;
//       int i = sorted[atomi].unsorted;
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
//       real uix = uind[i][0];
//       real uiy = uind[i][1];
//       real uiz = uind[i][2];
//       real corei = pcore[i];
//       real alphai = palpha[i];
//       real vali = pval[i];


//       int shatomk;
//       shatomk = lst[iw * WARP_SIZE + ilane];
//       real shx = sorted[shatomk].x;
//       real shy = sorted[shatomk].y;
//       real shz = sorted[shatomk].z;
//       int shk = sorted[shatomk].unsorted;
//       real shck = rpole[shk][mpl_pme_0];
//       real shdkx = rpole[shk][mpl_pme_x];
//       real shdky = rpole[shk][mpl_pme_y];
//       real shdkz = rpole[shk][mpl_pme_z];
//       real shqkxx = rpole[shk][mpl_pme_xx];
//       real shqkxy = rpole[shk][mpl_pme_xy];
//       real shqkxz = rpole[shk][mpl_pme_xz];
//       real shqkyy = rpole[shk][mpl_pme_yy];
//       real shqkyz = rpole[shk][mpl_pme_yz];
//       real shqkzz = rpole[shk][mpl_pme_zz];
//       real shukx = uind[shk][0];
//       real shuky = uind[shk][1];
//       real shukz = uind[shk][2];
//       real shcorek = pcore[shk];
//       real shalphak = palpha[shk];
//       real shvalk = pval[shk];


//       for (int j = 0; j < WARP_SIZE; ++j) {
//          int srclane = (ilane + j) & (WARP_SIZE - 1);
//          int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
//          real xr = __shfl_sync(ALL_LANES, shx, srclane) - xi;
//          real yr = __shfl_sync(ALL_LANES, shy, srclane) - yi;
//          real zr = __shfl_sync(ALL_LANES, shz, srclane) - zi;
//          int k = __shfl_sync(ALL_LANES, shk, srclane);
//          real ck = __shfl_sync(ALL_LANES, shck, srclane);
//          real dkx = __shfl_sync(ALL_LANES, shdkx, srclane);
//          real dky = __shfl_sync(ALL_LANES, shdky, srclane);
//          real dkz = __shfl_sync(ALL_LANES, shdkz, srclane);
//          real qkxx = __shfl_sync(ALL_LANES, shqkxx, srclane);
//          real qkxy = __shfl_sync(ALL_LANES, shqkxy, srclane);
//          real qkxz = __shfl_sync(ALL_LANES, shqkxz, srclane);
//          real qkyy = __shfl_sync(ALL_LANES, shqkyy, srclane);
//          real qkyz = __shfl_sync(ALL_LANES, shqkyz, srclane);
//          real qkzz = __shfl_sync(ALL_LANES, shqkzz, srclane);
//          real ukx = __shfl_sync(ALL_LANES, shukx, srclane);
//          real uky = __shfl_sync(ALL_LANES, shuky, srclane);
//          real ukz = __shfl_sync(ALL_LANES, shukz, srclane);
//          real corek = __shfl_sync(ALL_LANES, shcorek, srclane);
//          real alphak = __shfl_sync(ALL_LANES, shalphak, srclane);
//          real valk = __shfl_sync(ALL_LANES, shvalk, srclane);


//          real e;
//          PairPolarGrad pgrad;
//          zero(pgrad);

//          real r2 = image2(xr, yr, zr);
//          if (atomi < atomk && r2 <= off2) {
//             if CONSTEXPR (eq<ETYP, EWALD>()) {
//                pair_polar_chgpen<do_e, do_g, EWALD>( //
//                   r2, xr, yr, zr, 1, 1,       //
//                   ci, dix, diy, diz, corei, vali, alphai,
//                   qixx, qixy, qixz, qiyy, qiyz, qizz, uix,
//                   uiy, uiz,  //
//                   ck, dkx, dky, dkz, corek, valk, alphak,
//                   qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukx,
//                   uky, ukz,  //
//                   f, aewald, e, pgrad);

//                // printf("%5.2f %5.2f %14.8f\n", alphai, alphak, r2);

//             }
//             if CONSTEXPR (eq<ETYP, NON_EWALD>()) {
//                pair_polar_chgpen<do_e, do_g, NON_EWALD>( //
//                   r2, xr, yr, zr, 1, 1,           //
//                   ci, dix, diy, diz, corei, vali, alphai,
//                   qixx, qixy, qixz, qiyy, qiyz, qizz, uix,
//                   uiy, uiz,  //
//                   ck, dkx, dky, dkz, corek, valk, alphak,
//                   qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukx,
//                   uky, ukz,  //
//                   f, 0, e, pgrad);
//             }


//             if CONSTEXPR (do_a)
//                ctl += 1;
//             if CONSTEXPR (do_e)
//                etl += e;
//             if CONSTEXPR (do_v) {
//                vtlxx += -xr * pgrad.frcx;
//                vtlxy += -0.5f * (yr * pgrad.frcx + xr * pgrad.frcy);
//                vtlxz += -0.5f * (zr * pgrad.frcx + xr * pgrad.frcz);
//                vtlyy += -yr * pgrad.frcy;
//                vtlyz += -0.5f * (zr * pgrad.frcy + yr * pgrad.frcz);
//                vtlzz += -zr * pgrad.frcz;
//             }
//          } // end if (include)


//          if CONSTEXPR (do_g) {
//             gxi += pgrad.frcx;
//             gyi += pgrad.frcy;
//             gzi += pgrad.frcz;
//             gxk[srclane + (threadIdx.x - ilane)] -= pgrad.frcx;
//             gyk[srclane + (threadIdx.x - ilane)] -= pgrad.frcy;
//             gzk[srclane + (threadIdx.x - ilane)] -= pgrad.frcz;


//             txi += pgrad.ufldi[0];
//             tyi += pgrad.ufldi[1];
//             tzi += pgrad.ufldi[2];
//             txk[srclane + (threadIdx.x - ilane)] += pgrad.ufldk[0];
//             tyk[srclane + (threadIdx.x - ilane)] += pgrad.ufldk[1];
//             tzk[srclane + (threadIdx.x - ilane)] += pgrad.ufldk[2];


//             dui[0] += pgrad.dufldi[0];
//             dui[1] += pgrad.dufldi[1];
//             dui[2] += pgrad.dufldi[2];
//             dui[3] += pgrad.dufldi[3];
//             dui[4] += pgrad.dufldi[4];
//             dui[5] += pgrad.dufldi[5];
//             duk[srclane + (threadIdx.x - ilane)][0] += pgrad.dufldk[0];
//             duk[srclane + (threadIdx.x - ilane)][1] += pgrad.dufldk[1];
//             duk[srclane + (threadIdx.x - ilane)][2] += pgrad.dufldk[2];
//             duk[srclane + (threadIdx.x - ilane)][3] += pgrad.dufldk[3];
//             duk[srclane + (threadIdx.x - ilane)][4] += pgrad.dufldk[4];
//             duk[srclane + (threadIdx.x - ilane)][5] += pgrad.dufldk[5];
//          }
//       } // end for (j)


//       if CONSTEXPR (do_a)
//          atomic_add(ctl, nep, offset);
//       if CONSTEXPR (do_e)
//          atomic_add(etl, ep, offset);
//       if CONSTEXPR (do_g) {
//          atomic_add(gxi, gx, i);
//          atomic_add(gyi, gy, i);
//          atomic_add(gzi, gz, i);
//          atomic_add(txi, &ufld[i][0]);
//          atomic_add(tyi, &ufld[i][1]);
//          atomic_add(tzi, &ufld[i][2]);
//          atomic_add(gxk[threadIdx.x], gx, shk);
//          atomic_add(gyk[threadIdx.x], gy, shk);
//          atomic_add(gzk[threadIdx.x], gz, shk);
//          atomic_add(txk[threadIdx.x], &ufld[shk][0]);
//          atomic_add(tyk[threadIdx.x], &ufld[shk][1]);
//          atomic_add(tzk[threadIdx.x], &ufld[shk][2]);


//          atomic_add(dui[0], &dufld[i][0]);
//          atomic_add(dui[1], &dufld[i][1]);
//          atomic_add(dui[2], &dufld[i][2]);
//          atomic_add(dui[3], &dufld[i][3]);
//          atomic_add(dui[4], &dufld[i][4]);
//          atomic_add(dui[5], &dufld[i][5]);
//          atomic_add(duk[threadIdx.x][0], &dufld[shk][0]);
//          atomic_add(duk[threadIdx.x][1], &dufld[shk][1]);
//          atomic_add(duk[threadIdx.x][2], &dufld[shk][2]);
//          atomic_add(duk[threadIdx.x][3], &dufld[shk][3]);
//          atomic_add(duk[threadIdx.x][4], &dufld[shk][4]);
//          atomic_add(duk[threadIdx.x][5], &dufld[shk][5]);
//       }
//       if CONSTEXPR (do_v)
//          atomic_add(vtlxx, vtlxy, vtlxz, vtlyy, vtlyz, vtlzz, vir_ep,
//          offset);
//    } // end for (iw)
// }


// template <class Ver>
// __global__
// void epolar_chgpen_cu2(POLARPARAS, const real* restrict x,
//                        const real* restrict y, const real* restrict z,
//                        int nmdwexclude, const int (*restrict mdwexclude)[2],
//                        const real (*restrict mdwexclude_scale)[3])
// {
//    constexpr bool do_e = Ver::e;
//    constexpr bool do_a = Ver::a;
//    constexpr bool do_g = Ver::g;
//    constexpr bool do_v = Ver::v;


//    for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < nmdwexclude;
//         ii += blockDim.x * gridDim.x) {
//       int offset = ii & (bufsize - 1);


//       int i = mdwexclude[ii][0];
//       int k = mdwexclude[ii][1];

//       real dscale = mdwexclude_scale[ii][1];
//       real wscale = mdwexclude_scale[ii][2];


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
//       real uix = uind[i][0];
//       real uiy = uind[i][1];
//       real uiz = uind[i][2];
//       real corei = pcore[i];
//       real alphai = palpha[i];
//       real vali = pval[i];


//       real xr = x[k] - xi;
//       real yr = y[k] - yi;
//       real zr = z[k] - zi;


//       real r2 = image2(xr, yr, zr);
//       if (r2 <= off2) {
//          real ck = rpole[k][mpl_pme_0];
//          real dkx = rpole[k][mpl_pme_x];
//          real dky = rpole[k][mpl_pme_y];
//          real dkz = rpole[k][mpl_pme_z];
//          real qkxx = rpole[k][mpl_pme_xx];
//          real qkxy = rpole[k][mpl_pme_xy];
//          real qkxz = rpole[k][mpl_pme_xz];
//          real qkyy = rpole[k][mpl_pme_yy];
//          real qkyz = rpole[k][mpl_pme_yz];
//          real qkzz = rpole[k][mpl_pme_zz];
//          real ukx = uind[k][0];
//          real uky = uind[k][1];
//          real ukz = uind[k][2];
//          real corek = pcore[k];
//          real alphak = palpha[k];
//          real valk = pval[k];


//          real e;
//          PairPolarGrad pgrad;
//          pair_polar_chgpen<do_e, do_g, NON_EWALD>( //
//             r2, xr, yr, zr, dscale, wscale,        //
//             ci, dix, diy, diz, corei, vali, alphai,
//             qixx, qixy, qixz, qiyy, qiyz, qizz, uix, uiy,
//             uiz,  //
//             ck, dkx, dky, dkz, corek, valk, alphak,
//             qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukx, uky,
//             ukz,  //
//             f, 0, e, pgrad);

//          // printf("%5.2f %5.2f %14.8f\n", alphai, alphak, r2);

//          if CONSTEXPR (do_a)
//             if (dscale == -1)
//                atomic_add(-1, nep, offset);
//          if CONSTEXPR (do_e)
//             atomic_add(e, ep, offset);
//          if CONSTEXPR (do_g) {
//             atomic_add(pgrad.frcx, gx, i);
//             atomic_add(pgrad.frcy, gy, i);
//             atomic_add(pgrad.frcz, gz, i);
//             atomic_add(-pgrad.frcx, gx, k);
//             atomic_add(-pgrad.frcy, gy, k);
//             atomic_add(-pgrad.frcz, gz, k);


//             atomic_add(pgrad.ufldi[0], &ufld[i][0]);
//             atomic_add(pgrad.ufldi[1], &ufld[i][1]);
//             atomic_add(pgrad.ufldi[2], &ufld[i][2]);
//             atomic_add(pgrad.ufldk[0], &ufld[k][0]);
//             atomic_add(pgrad.ufldk[1], &ufld[k][1]);
//             atomic_add(pgrad.ufldk[2], &ufld[k][2]);


//             atomic_add(pgrad.dufldi[0], &dufld[i][0]);
//             atomic_add(pgrad.dufldi[1], &dufld[i][1]);
//             atomic_add(pgrad.dufldi[2], &dufld[i][2]);
//             atomic_add(pgrad.dufldi[3], &dufld[i][3]);
//             atomic_add(pgrad.dufldi[4], &dufld[i][4]);
//             atomic_add(pgrad.dufldi[5], &dufld[i][5]);
//             atomic_add(pgrad.dufldk[0], &dufld[k][0]);
//             atomic_add(pgrad.dufldk[1], &dufld[k][1]);
//             atomic_add(pgrad.dufldk[2], &dufld[k][2]);
//             atomic_add(pgrad.dufldk[3], &dufld[k][3]);
//             atomic_add(pgrad.dufldk[4], &dufld[k][4]);
//             atomic_add(pgrad.dufldk[5], &dufld[k][5]);


//             if CONSTEXPR (do_v) {
//                real vxx = -xr * pgrad.frcx;
//                real vxy = -0.5f * (yr * pgrad.frcx + xr * pgrad.frcy);
//                real vxz = -0.5f * (zr * pgrad.frcx + xr * pgrad.frcz);
//                real vyy = -yr * pgrad.frcy;
//                real vyz = -0.5f * (zr * pgrad.frcy + yr * pgrad.frcz);
//                real vzz = -zr * pgrad.frcz;
//                atomic_add(vxx, vxy, vxz, vyy, vyz, vzz, vir_ep, offset);
//             }
//          }
//       } // end if (r2 <= off2)
//    }
// }

template <class Ver, class ETYP>
__global__
void epolar_chgpen_cu1(
   int n, TINKER_IMAGE_PARAMS, count_buffer restrict np,
   energy_buffer restrict ep, virial_buffer restrict vp, grad_prec* restrict gx,
   grad_prec* restrict gy, grad_prec* restrict gz, real off,
   const unsigned* restrict dwinfo, int nexclude,
   const int (*restrict exclude)[2], const real (*restrict exclude_scale)[3],
   const real* restrict x, const real* restrict y, const real* restrict z,
   const Spatial::SortedAtom* restrict sorted, int nakpl,
   const int* restrict iakpl, int niak, const int* restrict iak,
   const int* restrict lst, real (*restrict ufld)[3], real (*restrict dufld)[6],
   const real (*restrict uind)[3], const real (*restrict rpole)[10],
   real* restrict pcore, real* restrict pval, const real* restrict palpha,
   real aewald, real f)
{
   constexpr bool do_a = Ver::a;
   constexpr bool do_e = Ver::e;
   constexpr bool do_v = Ver::v;
   constexpr bool do_g = Ver::g;


   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);


   int nptl;
   if CONSTEXPR (do_a) {
      nptl = 0;
   }
   using ebuf_prec = energy_buffer_traits::type;
   ebuf_prec eptl;
   if CONSTEXPR (do_e) {
      eptl = 0;
   }
   using vbuf_prec = virial_buffer_traits::type;
   vbuf_prec vptlxx, vptlyx, vptlzx, vptlyy, vptlzy, vptlzz;
   if CONSTEXPR (do_v) {
      vptlxx = 0;
      vptlyx = 0;
      vptlzx = 0;
      vptlyy = 0;
      vptlzy = 0;
      vptlzz = 0;
   }


   __shared__ real shxi[BLOCK_DIM];
   __shared__ real shyi[BLOCK_DIM];
   __shared__ real shzi[BLOCK_DIM];
   real xk;
   real yk;
   real zk;
   __shared__ real shgxi[BLOCK_DIM];
   __shared__ real shgyi[BLOCK_DIM];
   __shared__ real shgzi[BLOCK_DIM];
   __shared__ real shtxi[BLOCK_DIM];
   __shared__ real shtyi[BLOCK_DIM];
   __shared__ real shtzi[BLOCK_DIM];
   __shared__ real shdui0[BLOCK_DIM];
   __shared__ real shdui1[BLOCK_DIM];
   __shared__ real shdui2[BLOCK_DIM];
   __shared__ real shdui3[BLOCK_DIM];
   __shared__ real shdui4[BLOCK_DIM];
   __shared__ real shdui5[BLOCK_DIM];
   real gxk;
   real gyk;
   real gzk;
   real txk;
   real tyk;
   real tzk;
   real duk0;
   real duk1;
   real duk2;
   real duk3;
   real duk4;
   real duk5;
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
   __shared__ real shuix[BLOCK_DIM];
   __shared__ real shuiy[BLOCK_DIM];
   __shared__ real shuiz[BLOCK_DIM];
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
   real ukx;
   real uky;
   real ukz;
   real corek;
   real alphak;
   real valk;


   //* /
   // exclude
   for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
      if CONSTEXPR (do_g) {
         shgxi[threadIdx.x] = 0;
         shgyi[threadIdx.x] = 0;
         shgzi[threadIdx.x] = 0;
         shtxi[threadIdx.x] = 0;
         shtyi[threadIdx.x] = 0;
         shtzi[threadIdx.x] = 0;
         shdui0[threadIdx.x] = 0;
         shdui1[threadIdx.x] = 0;
         shdui2[threadIdx.x] = 0;
         shdui3[threadIdx.x] = 0;
         shdui4[threadIdx.x] = 0;
         shdui5[threadIdx.x] = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
         txk = 0;
         tyk = 0;
         tzk = 0;
         duk0 = 0;
         duk1 = 0;
         duk2 = 0;
         duk3 = 0;
         duk4 = 0;
         duk5 = 0;
      }


      int shi = exclude[ii][0];
      int k = exclude[ii][1];
      real scaleb = exclude_scale[ii][1];
      real scalec = exclude_scale[ii][2];


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
      real uix = uind[shi][0];
      real uiy = uind[shi][1];
      real uiz = uind[shi][2];
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
      ukx = uind[k][0];
      uky = uind[k][1];
      ukz = uind[k][2];
      corek = pcore[k];
      alphak = palpha[k];
      valk = pval[k];


      int klane = threadIdx.x;
      constexpr bool incl = true;
      real xr = xk - xi;
      real yr = yk - yi;
      real zr = zk - zi;

      real e;
      PairPolarGrad pgrad;
      zero(pgrad);

      real r2 = image2(xr, yr, zr);
      if (r2 <= off * off and incl) {
         pair_polar_chgpen<do_e, do_g, ETYP>( //
            r2, xr, yr, zr, scaleb, scalec,   //
            ci, dix, diy, diz, corei, vali, alphai, qixx, qixy, qixz, qiyy,
            qiyz, qizz, uix, uiy, uiz, //
            ck, dkx, dky, dkz, corek, valk, alphak, qkxx, qkxy, qkxz, qkyy,
            qkyz, qkzz, ukx, uky, ukz, //
            f, aewald, e, pgrad);

            //printf("block1 %2d %2d %5.2f %5.2f %5.2f\n", shi+1, k+1, yi, yk, yr);

         if CONSTEXPR (do_a)
            if (e != 0 and (scaleb + scalec) != 0)
               nptl += 1;
         if CONSTEXPR (do_e)
            eptl += e;
         if CONSTEXPR (do_g) {
            shgxi[klane] += pgrad.frcx;
            shgyi[klane] += pgrad.frcy;
            shgzi[klane] += pgrad.frcz;
            gxk -= pgrad.frcx;
            gyk -= pgrad.frcy;
            gzk -= pgrad.frcz;


            shtxi[klane] += pgrad.ufldi[0];
            shtyi[klane] += pgrad.ufldi[1];
            shtzi[klane] += pgrad.ufldi[2];
            txk += pgrad.ufldk[0];
            tyk += pgrad.ufldk[1];
            tzk += pgrad.ufldk[2];


            shdui0[klane] += pgrad.dufldi[0];
            shdui1[klane] += pgrad.dufldi[1];
            shdui2[klane] += pgrad.dufldi[2];
            shdui3[klane] += pgrad.dufldi[3];
            shdui4[klane] += pgrad.dufldi[4];
            shdui5[klane] += pgrad.dufldi[5];
            duk0 += pgrad.dufldk[0];
            duk1 += pgrad.dufldk[1];
            duk2 += pgrad.dufldk[2];
            duk3 += pgrad.dufldk[3];
            duk4 += pgrad.dufldk[4];
            duk5 += pgrad.dufldk[5];
         }
         if CONSTEXPR (do_v) {
            vptlxx += -xr * pgrad.frcx;
            vptlyx += -0.5f * (yr * pgrad.frcx + xr * pgrad.frcy);
            vptlzx += -0.5f * (zr * pgrad.frcx + xr * pgrad.frcz);
            vptlyy += -yr * pgrad.frcy;
            vptlzy += -0.5f * (zr * pgrad.frcy + yr * pgrad.frcz);
            vptlzz += -zr * pgrad.frcz;
         }
      } // end if (include)


      if CONSTEXPR (do_g) {
         atomic_add(shgxi[threadIdx.x], gx, shi);
         atomic_add(shgyi[threadIdx.x], gy, shi);
         atomic_add(shgzi[threadIdx.x], gz, shi);
         atomic_add(shtxi[threadIdx.x], &ufld[shi][0]);
         atomic_add(shtyi[threadIdx.x], &ufld[shi][1]);
         atomic_add(shtzi[threadIdx.x], &ufld[shi][2]);
         atomic_add(shdui0[threadIdx.x], &dufld[shi][0]);
         atomic_add(shdui1[threadIdx.x], &dufld[shi][1]);
         atomic_add(shdui2[threadIdx.x], &dufld[shi][2]);
         atomic_add(shdui3[threadIdx.x], &dufld[shi][3]);
         atomic_add(shdui4[threadIdx.x], &dufld[shi][4]);
         atomic_add(shdui5[threadIdx.x], &dufld[shi][5]);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
         atomic_add(txk, &ufld[k][0]);
         atomic_add(tyk, &ufld[k][1]);
         atomic_add(tzk, &ufld[k][2]);
         atomic_add(duk0, &dufld[k][0]);
         atomic_add(duk1, &dufld[k][1]);
         atomic_add(duk2, &dufld[k][2]);
         atomic_add(duk3, &dufld[k][3]);
         atomic_add(duk4, &dufld[k][4]);
         atomic_add(duk5, &dufld[k][5]);
      }
   }
   // */


   //* /
   // block pairs that have scale factors
   for (int iw = iwarp; iw < nakpl; iw += nwarp) {
      if CONSTEXPR (do_g) {
         shgxi[threadIdx.x] = 0;
         shgyi[threadIdx.x] = 0;
         shgzi[threadIdx.x] = 0;
         shtxi[threadIdx.x] = 0;
         shtyi[threadIdx.x] = 0;
         shtzi[threadIdx.x] = 0;
         shdui0[threadIdx.x] = 0;
         shdui1[threadIdx.x] = 0;
         shdui2[threadIdx.x] = 0;
         shdui3[threadIdx.x] = 0;
         shdui4[threadIdx.x] = 0;
         shdui5[threadIdx.x] = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
         txk = 0;
         tyk = 0;
         tzk = 0;
         duk0 = 0;
         duk1 = 0;
         duk2 = 0;
         duk3 = 0;
         duk4 = 0;
         duk5 = 0;
      }


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
      shuix[threadIdx.x] = uind[shi][0];
      shuiy[threadIdx.x] = uind[shi][1];
      shuiz[threadIdx.x] = uind[shi][2];
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
      ukx = uind[k][0];
      uky = uind[k][1];
      ukz = uind[k][2];
      corek = pcore[k];
      alphak = palpha[k];
      valk = pval[k];


      unsigned int dwinfo0 = dwinfo[iw * WARP_SIZE + ilane];


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
         real uix = shuix[klane];
         real uiy = shuiy[klane];
         real uiz = shuiz[klane];
         real corei = shcorei[klane];
         real alphai = shalphai[klane];
         real vali = shvali[klane];


         bool incl = iid < kid and kid < n;
         incl = incl and (dwinfo0 & srcmask) == 0;
         real scaleb = 1;
         real scalec = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;

         real e;
         PairPolarGrad pgrad;
         zero(pgrad);

         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_polar_chgpen<do_e, do_g, ETYP>( //
               r2, xr, yr, zr, scaleb, scalec,   //
               ci, dix, diy, diz, corei, vali, alphai, qixx, qixy, qixz, qiyy,
               qiyz, qizz, uix, uiy, uiz, //
               ck, dkx, dky, dkz, corek, valk, alphak, qkxx, qkxy, qkxz, qkyy,
               qkyz, qkzz, ukx, uky, ukz, //
               f, aewald, e, pgrad);

            //printf("block2 %2d %2d %5.2f %5.2f %5.2f\n", iid+1, kid+1, yi, yk, yr);

            if CONSTEXPR (do_a)
               if (e != 0)
                  nptl += 1;
            if CONSTEXPR (do_e)
               eptl += e;
            if CONSTEXPR (do_g) {
               shgxi[klane] += pgrad.frcx;
               shgyi[klane] += pgrad.frcy;
               shgzi[klane] += pgrad.frcz;
               gxk -= pgrad.frcx;
               gyk -= pgrad.frcy;
               gzk -= pgrad.frcz;


               shtxi[klane] += pgrad.ufldi[0];
               shtyi[klane] += pgrad.ufldi[1];
               shtzi[klane] += pgrad.ufldi[2];
               txk += pgrad.ufldk[0];
               tyk += pgrad.ufldk[1];
               tzk += pgrad.ufldk[2];


               shdui0[klane] += pgrad.dufldi[0];
               shdui1[klane] += pgrad.dufldi[1];
               shdui2[klane] += pgrad.dufldi[2];
               shdui3[klane] += pgrad.dufldi[3];
               shdui4[klane] += pgrad.dufldi[4];
               shdui5[klane] += pgrad.dufldi[5];
               duk0 += pgrad.dufldk[0];
               duk1 += pgrad.dufldk[1];
               duk2 += pgrad.dufldk[2];
               duk3 += pgrad.dufldk[3];
               duk4 += pgrad.dufldk[4];
               duk5 += pgrad.dufldk[5];
            }
            if CONSTEXPR (do_v) {
               vptlxx += -xr * pgrad.frcx;
               vptlyx += -0.5f * (yr * pgrad.frcx + xr * pgrad.frcy);
               vptlzx += -0.5f * (zr * pgrad.frcx + xr * pgrad.frcz);
               vptlyy += -yr * pgrad.frcy;
               vptlzy += -0.5f * (zr * pgrad.frcy + yr * pgrad.frcz);
               vptlzz += -zr * pgrad.frcz;
            }
         } // end if (include)


         shiid = __shfl_sync(ALL_LANES, shiid, ilane + 1);
      }


      if CONSTEXPR (do_g) {
         atomic_add(shgxi[threadIdx.x], gx, shi);
         atomic_add(shgyi[threadIdx.x], gy, shi);
         atomic_add(shgzi[threadIdx.x], gz, shi);
         atomic_add(shtxi[threadIdx.x], &ufld[shi][0]);
         atomic_add(shtyi[threadIdx.x], &ufld[shi][1]);
         atomic_add(shtzi[threadIdx.x], &ufld[shi][2]);
         atomic_add(shdui0[threadIdx.x], &dufld[shi][0]);
         atomic_add(shdui1[threadIdx.x], &dufld[shi][1]);
         atomic_add(shdui2[threadIdx.x], &dufld[shi][2]);
         atomic_add(shdui3[threadIdx.x], &dufld[shi][3]);
         atomic_add(shdui4[threadIdx.x], &dufld[shi][4]);
         atomic_add(shdui5[threadIdx.x], &dufld[shi][5]);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
         atomic_add(txk, &ufld[k][0]);
         atomic_add(tyk, &ufld[k][1]);
         atomic_add(tzk, &ufld[k][2]);
         atomic_add(duk0, &dufld[k][0]);
         atomic_add(duk1, &dufld[k][1]);
         atomic_add(duk2, &dufld[k][2]);
         atomic_add(duk3, &dufld[k][3]);
         atomic_add(duk4, &dufld[k][4]);
         atomic_add(duk5, &dufld[k][5]);
      }
   }
   // */


   //* /
   // block-atoms
   for (int iw = iwarp; iw < niak; iw += nwarp) {
      if CONSTEXPR (do_g) {
         shgxi[threadIdx.x] = 0;
         shgyi[threadIdx.x] = 0;
         shgzi[threadIdx.x] = 0;
         shtxi[threadIdx.x] = 0;
         shtyi[threadIdx.x] = 0;
         shtzi[threadIdx.x] = 0;
         shdui0[threadIdx.x] = 0;
         shdui1[threadIdx.x] = 0;
         shdui2[threadIdx.x] = 0;
         shdui3[threadIdx.x] = 0;
         shdui4[threadIdx.x] = 0;
         shdui5[threadIdx.x] = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
         txk = 0;
         tyk = 0;
         tzk = 0;
         duk0 = 0;
         duk1 = 0;
         duk2 = 0;
         duk3 = 0;
         duk4 = 0;
         duk5 = 0;
      }


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
      shuix[threadIdx.x] = uind[shi][0];
      shuiy[threadIdx.x] = uind[shi][1];
      shuiz[threadIdx.x] = uind[shi][2];
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
      ukx = uind[k][0];
      uky = uind[k][1];
      ukz = uind[k][2];
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
         real uix = shuix[klane];
         real uiy = shuiy[klane];
         real uiz = shuiz[klane];
         real corei = shcorei[klane];
         real alphai = shalphai[klane];
         real vali = shvali[klane];


         bool incl = atomk > 0;
         real scaleb = 1;
         real scalec = 1;
         real xr = xk - xi;
         real yr = yk - yi;
         real zr = zk - zi;

         real e;
         PairPolarGrad pgrad;
         zero(pgrad);

         real r2 = image2(xr, yr, zr);
         if (r2 <= off * off and incl) {
            pair_polar_chgpen<do_e, do_g, ETYP>( //
               r2, xr, yr, zr, scaleb, scalec,   //
               ci, dix, diy, diz, corei, vali, alphai, qixx, qixy, qixz, qiyy,
               qiyz, qizz, uix, uiy, uiz, //
               ck, dkx, dky, dkz, corek, valk, alphak, qkxx, qkxy, qkxz, qkyy,
               qkyz, qkzz, ukx, uky, ukz, //
               f, aewald, e, pgrad);


            if CONSTEXPR (do_a)
               if (e != 0)
                  nptl += 1;
            if CONSTEXPR (do_e)
               eptl += e;
            if CONSTEXPR (do_g) {
               shgxi[klane] += pgrad.frcx;
               shgyi[klane] += pgrad.frcy;
               shgzi[klane] += pgrad.frcz;
               gxk -= pgrad.frcx;
               gyk -= pgrad.frcy;
               gzk -= pgrad.frcz;


               shtxi[klane] += pgrad.ufldi[0];
               shtyi[klane] += pgrad.ufldi[1];
               shtzi[klane] += pgrad.ufldi[2];
               txk += pgrad.ufldk[0];
               tyk += pgrad.ufldk[1];
               tzk += pgrad.ufldk[2];


               shdui0[klane] += pgrad.dufldi[0];
               shdui1[klane] += pgrad.dufldi[1];
               shdui2[klane] += pgrad.dufldi[2];
               shdui3[klane] += pgrad.dufldi[3];
               shdui4[klane] += pgrad.dufldi[4];
               shdui5[klane] += pgrad.dufldi[5];
               duk0 += pgrad.dufldk[0];
               duk1 += pgrad.dufldk[1];
               duk2 += pgrad.dufldk[2];
               duk3 += pgrad.dufldk[3];
               duk4 += pgrad.dufldk[4];
               duk5 += pgrad.dufldk[5];
            }
            if CONSTEXPR (do_v) {
               vptlxx += -xr * pgrad.frcx;
               vptlyx += -0.5f * (yr * pgrad.frcx + xr * pgrad.frcy);
               vptlzx += -0.5f * (zr * pgrad.frcx + xr * pgrad.frcz);
               vptlyy += -yr * pgrad.frcy;
               vptlzy += -0.5f * (zr * pgrad.frcy + yr * pgrad.frcz);
               vptlzz += -zr * pgrad.frcz;
            }
         } // end if (include)
      }


      if CONSTEXPR (do_g) {
         atomic_add(shgxi[threadIdx.x], gx, shi);
         atomic_add(shgyi[threadIdx.x], gy, shi);
         atomic_add(shgzi[threadIdx.x], gz, shi);
         atomic_add(shtxi[threadIdx.x], &ufld[shi][0]);
         atomic_add(shtyi[threadIdx.x], &ufld[shi][1]);
         atomic_add(shtzi[threadIdx.x], &ufld[shi][2]);
         atomic_add(shdui0[threadIdx.x], &dufld[shi][0]);
         atomic_add(shdui1[threadIdx.x], &dufld[shi][1]);
         atomic_add(shdui2[threadIdx.x], &dufld[shi][2]);
         atomic_add(shdui3[threadIdx.x], &dufld[shi][3]);
         atomic_add(shdui4[threadIdx.x], &dufld[shi][4]);
         atomic_add(shdui5[threadIdx.x], &dufld[shi][5]);
         atomic_add(gxk, gx, k);
         atomic_add(gyk, gy, k);
         atomic_add(gzk, gz, k);
         atomic_add(txk, &ufld[k][0]);
         atomic_add(tyk, &ufld[k][1]);
         atomic_add(tzk, &ufld[k][2]);
         atomic_add(duk0, &dufld[k][0]);
         atomic_add(duk1, &dufld[k][1]);
         atomic_add(duk2, &dufld[k][2]);
         atomic_add(duk3, &dufld[k][3]);
         atomic_add(duk4, &dufld[k][4]);
         atomic_add(duk5, &dufld[k][5]);
      }
   }
   // */


   if CONSTEXPR (do_a) {
      atomic_add(nptl, np, ithread);
   }
   if CONSTEXPR (do_e) {
      atomic_add(eptl, ep, ithread);
   }
   if CONSTEXPR (do_v) {
      atomic_add(vptlxx, vptlyx, vptlzx, vptlyy, vptlzy, vptlzz, vp, ithread);
   }
} // generated by ComplexKernelBuilder (ck.py) 1.5.1


template <class Ver, class ETYP>
void epolar_chgpen_cu(const real (*uind)[3])
{
   constexpr bool do_g = Ver::g;


   // const auto& st = *mspatial_unit;
   const auto& st = *mspatial_v2_unit;
   real off;

   if CONSTEXPR (eq<ETYP, EWALD>())
      off = switch_off(switch_ewald);
   else
      off = switch_off(switch_mpole);
   auto bufsize = buffer_size();


   const real f = 0.5f * electric / dielec;
   real aewald = 0;
   if CONSTEXPR (eq<ETYP, EWALD>()) {
      PMEUnit pu = ppme_unit;
      aewald = pu->aewald;
   }


   if CONSTEXPR (do_g) {
      darray::zero(PROCEED_NEW_Q, n, ufld, dufld);
   }
   // if (st.niak > 0) {
   //    auto ker1 = epolar_chgpen_cu1<Ver, ETYP>;
   //    launch_k1s(nonblk, WARP_SIZE * st.niak, ker1, //
   //               bufsize, nep, ep, vir_ep, depx, depy, depz, ufld, dufld,
   //               TINKER_IMAGE_ARGS, off2, f, rpole, pcore, pval, palpha,
   //               uind, st.sorted, st.niak, st.iak, st.lst, n, aewald);
   // }
   // if (nmdwexclude > 0) {
   //    auto ker2 = epolar_chgpen_cu2<Ver>;
   //    launch_k1s(nonblk, nmdwexclude, ker2, //
   //               bufsize, nep, ep, vir_ep, depx, depy, depz, ufld, dufld,
   //               TINKER_IMAGE_ARGS, off2, f, rpole, pcore, pval, palpha,
   //               uind, //
   //               x, y, z, nmdwexclude, mdwexclude, mdwexclude_scale);
   // }

   int ngrid = get_grid_size(BLOCK_DIM);
   epolar_chgpen_cu1<Ver, ETYP><<<ngrid, BLOCK_DIM, 0, nonblk>>>(
      st.n, TINKER_IMAGE_ARGS, nep, ep, vir_ep, depx, depy, depz, off,
      st.si1.bit0, nmdwexclude, mdwexclude, mdwexclude_scale, st.x, st.y, st.z,
      st.sorted, st.nakpl, st.iakpl, st.niak, st.iak, st.lst, ufld, dufld, uind,
      rpole, pcore, pval, palpha, aewald, f);

   // torque
   if CONSTEXPR (do_g) {
      launch_k1s(nonblk, n, epolar_trq_cu, //
                 trqx, trqy, trqz, n, rpole, ufld, dufld);
   }
}


void epolar_chgpen_nonewald_cu(int vers, const real (*uind)[3])
{
   if (vers == calc::v0) {
      epolar_chgpen_cu<calc::V0, NON_EWALD>(uind);
   } else if (vers == calc::v1) {
      epolar_chgpen_cu<calc::V1, NON_EWALD>(uind);
   } else if (vers == calc::v3) {
      epolar_chgpen_cu<calc::V3, NON_EWALD>(uind);
   } else if (vers == calc::v4) {
      epolar_chgpen_cu<calc::V4, NON_EWALD>(uind);
   } else if (vers == calc::v5) {
      epolar_chgpen_cu<calc::V5, NON_EWALD>(uind);
   } else if (vers == calc::v6) {
      epolar_chgpen_cu<calc::V6, NON_EWALD>(uind);
   }
}


void epolar_chgpen_ewald_real_cu(int vers, const real (*uind)[3])
{
   if (vers == calc::v0) {
      epolar_chgpen_cu<calc::V0, EWALD>(uind);
   } else if (vers == calc::v1) {
      epolar_chgpen_cu<calc::V1, EWALD>(uind);
   } else if (vers == calc::v3) {
      epolar_chgpen_cu<calc::V3, EWALD>(uind);
   } else if (vers == calc::v4) {
      epolar_chgpen_cu<calc::V4, EWALD>(uind);
   } else if (vers == calc::v5) {
      epolar_chgpen_cu<calc::V5, EWALD>(uind);
   } else if (vers == calc::v6) {
      epolar_chgpen_cu<calc::V6, EWALD>(uind);
   }
}
}
