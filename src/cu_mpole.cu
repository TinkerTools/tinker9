#include "add.h"
#include "e_mpole.h"
#include "empole_self.h"
#include "launch.h"
#include "md.h"
#include "named_struct.h"
#include "pme.h"
#include "seq_image.h"
#include "seq_pair_mpole.h"
#include "spatial.h"


TINKER_NAMESPACE_BEGIN
#define EMPOLE_ARGS                                                            \
   size_t bufsize, count_buffer restrict nem, energy_buffer restrict em,       \
      virial_buffer restrict vir_em, grad_prec *restrict gx,                   \
      grad_prec *restrict gy, grad_prec *restrict gz, real *restrict trqx,     \
      real *restrict trqy, real *restrict trqz, TINKER_IMAGE_PARAMS,           \
      real off2, real f, const real(*restrict rpole)[10]


template <class Ver, class ETYP>
__global__
void empole_cu1(EMPOLE_ARGS, const Spatial::SortedAtom* restrict sorted,
                int niak, const int* restrict iak, const int* restrict lst,
                int n, real aewald)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);
   const int offset = ithread & (bufsize - 1);


   MAYBE_UNUSED int ctl;
   MAYBE_UNUSED real etl;
   MAYBE_UNUSED real vtlxx, vtlxy, vtlxz, vtlyy, vtlyz, vtlzz;
   MAYBE_UNUSED real gxi, gyi, gzi, txi, tyi, tzi;
   MAYBE_UNUSED __shared__ real gxk[BLOCK_DIM], gyk[BLOCK_DIM], gzk[BLOCK_DIM],
      txk[BLOCK_DIM], tyk[BLOCK_DIM], tzk[BLOCK_DIM];


   for (int iw = iwarp; iw < niak; iw += nwarp) {
      if CONSTEXPR (do_a)
         ctl = 0;
      if CONSTEXPR (do_e)
         etl = 0;
      if CONSTEXPR (do_v) {
         vtlxx = 0;
         vtlxy = 0;
         vtlxz = 0;
         vtlyy = 0;
         vtlyz = 0;
         vtlzz = 0;
      }
      if CONSTEXPR (do_g) {
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
      }


      int atomi;
      atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
      real xi = sorted[atomi].x;
      real yi = sorted[atomi].y;
      real zi = sorted[atomi].z;
      int i = sorted[atomi].unsorted;
      real ci = rpole[i][mpl_pme_0];
      real dix = rpole[i][mpl_pme_x];
      real diy = rpole[i][mpl_pme_y];
      real diz = rpole[i][mpl_pme_z];
      real qixx = rpole[i][mpl_pme_xx];
      real qixy = rpole[i][mpl_pme_xy];
      real qixz = rpole[i][mpl_pme_xz];
      real qiyy = rpole[i][mpl_pme_yy];
      real qiyz = rpole[i][mpl_pme_yz];
      real qizz = rpole[i][mpl_pme_zz];


      int shatomk;
      shatomk = lst[iw * WARP_SIZE + ilane];
      real shx = sorted[shatomk].x;
      real shy = sorted[shatomk].y;
      real shz = sorted[shatomk].z;
      int shk = sorted[shatomk].unsorted;
      real shck = rpole[shk][mpl_pme_0];
      real shdkx = rpole[shk][mpl_pme_x];
      real shdky = rpole[shk][mpl_pme_y];
      real shdkz = rpole[shk][mpl_pme_z];
      real shqkxx = rpole[shk][mpl_pme_xx];
      real shqkxy = rpole[shk][mpl_pme_xy];
      real shqkxz = rpole[shk][mpl_pme_xz];
      real shqkyy = rpole[shk][mpl_pme_yy];
      real shqkyz = rpole[shk][mpl_pme_yz];
      real shqkzz = rpole[shk][mpl_pme_zz];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
         real xr = __shfl_sync(ALL_LANES, shx, srclane) - xi;
         real yr = __shfl_sync(ALL_LANES, shy, srclane) - yi;
         real zr = __shfl_sync(ALL_LANES, shz, srclane) - zi;
         int k = __shfl_sync(ALL_LANES, shk, srclane);
         real ck = __shfl_sync(ALL_LANES, shck, srclane);
         real dkx = __shfl_sync(ALL_LANES, shdkx, srclane);
         real dky = __shfl_sync(ALL_LANES, shdky, srclane);
         real dkz = __shfl_sync(ALL_LANES, shdkz, srclane);
         real qkxx = __shfl_sync(ALL_LANES, shqkxx, srclane);
         real qkxy = __shfl_sync(ALL_LANES, shqkxy, srclane);
         real qkxz = __shfl_sync(ALL_LANES, shqkxz, srclane);
         real qkyy = __shfl_sync(ALL_LANES, shqkyy, srclane);
         real qkyz = __shfl_sync(ALL_LANES, shqkyz, srclane);
         real qkzz = __shfl_sync(ALL_LANES, shqkzz, srclane);


         PairMPoleGrad pgrad;
         zero(pgrad);


         real r2 = image2(xr, yr, zr);
         if (atomi < atomk && r2 <= off2) {
            real e = 0;
            if CONSTEXPR (eq<ETYP, EWALD>()) {
               pair_mpole<do_e, do_g, EWALD>(
                  r2, xr, yr, zr, 1,                                     //
                  ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz, //
                  ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, //
                  f, aewald, e, pgrad);
            }
            if CONSTEXPR (eq<ETYP, NON_EWALD>()) {
               pair_mpole<do_e, do_g, NON_EWALD>(
                  r2, xr, yr, zr, 1,                                     //
                  ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz, //
                  ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, //
                  f, 0, e, pgrad);
            }


            if CONSTEXPR (do_a)
               ctl += 1;
            if CONSTEXPR (do_e)
               etl += e;
            if CONSTEXPR (do_v) {
               vtlxx += -xr * pgrad.frcx;
               vtlxy += -0.5f * (yr * pgrad.frcx + xr * pgrad.frcy);
               vtlxz += -0.5f * (zr * pgrad.frcx + xr * pgrad.frcz);
               vtlyy += -yr * pgrad.frcy;
               vtlyz += -0.5f * (zr * pgrad.frcy + yr * pgrad.frcz);
               vtlzz += -zr * pgrad.frcz;
            }
         } // enf if (include)


         if CONSTEXPR (do_g) {
            gxi += pgrad.frcx;
            gyi += pgrad.frcy;
            gzi += pgrad.frcz;
            gxk[srclane + (threadIdx.x - ilane)] -= pgrad.frcx;
            gyk[srclane + (threadIdx.x - ilane)] -= pgrad.frcy;
            gzk[srclane + (threadIdx.x - ilane)] -= pgrad.frcz;


            txi += pgrad.ttmi[0];
            tyi += pgrad.ttmi[1];
            tzi += pgrad.ttmi[2];
            txk[srclane + (threadIdx.x - ilane)] += pgrad.ttmk[0];
            tyk[srclane + (threadIdx.x - ilane)] += pgrad.ttmk[1];
            tzk[srclane + (threadIdx.x - ilane)] += pgrad.ttmk[2];
         }
      } // end for (j)


      if CONSTEXPR (do_a)
         atomic_add(ctl, nem, offset);
      if CONSTEXPR (do_e)
         atomic_add(etl, em, offset);
      if CONSTEXPR (do_g) {
         atomic_add(gxi, gx, i);
         atomic_add(gyi, gy, i);
         atomic_add(gzi, gz, i);
         atomic_add(txi, trqx, i);
         atomic_add(tyi, trqy, i);
         atomic_add(tzi, trqz, i);
         atomic_add(gxk[threadIdx.x], gx, shk);
         atomic_add(gyk[threadIdx.x], gy, shk);
         atomic_add(gzk[threadIdx.x], gz, shk);
         atomic_add(txk[threadIdx.x], trqx, shk);
         atomic_add(tyk[threadIdx.x], trqy, shk);
         atomic_add(tzk[threadIdx.x], trqz, shk);
      }
      if CONSTEXPR (do_v)
         atomic_add(vtlxx, vtlxy, vtlxz, vtlyy, vtlyz, vtlzz, vir_em, offset);
   } // end for (iw)
}


template <class Ver>
__global__
void empole_cu2(EMPOLE_ARGS, const real* restrict x, const real* restrict y,
                const real* restrict z, int nmexclude_,
                const int (*restrict mexclude_)[2],
                const real* restrict mexclude_scale_)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < nmexclude_;
        ii += blockDim.x * gridDim.x) {
      int offset = ii & (bufsize - 1);


      int i = mexclude_[ii][0];
      int k = mexclude_[ii][1];
      real mscale = mexclude_scale_[ii];


      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real ci = rpole[i][mpl_pme_0];
      real dix = rpole[i][mpl_pme_x];
      real diy = rpole[i][mpl_pme_y];
      real diz = rpole[i][mpl_pme_z];
      real qixx = rpole[i][mpl_pme_xx];
      real qixy = rpole[i][mpl_pme_xy];
      real qixz = rpole[i][mpl_pme_xz];
      real qiyy = rpole[i][mpl_pme_yy];
      real qiyz = rpole[i][mpl_pme_yz];
      real qizz = rpole[i][mpl_pme_zz];


      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off2) {
         real e;
         PairMPoleGrad pgrad;
         pair_mpole<do_e, do_g, NON_EWALD>(
            r2, xr, yr, zr, mscale, ci, dix, diy, diz, qixx, qixy, qixz, qiyy,
            qiyz, qizz, rpole[k][mpl_pme_0], rpole[k][mpl_pme_x],
            rpole[k][mpl_pme_y], rpole[k][mpl_pme_z], rpole[k][mpl_pme_xx],
            rpole[k][mpl_pme_xy], rpole[k][mpl_pme_xz], rpole[k][mpl_pme_yy],
            rpole[k][mpl_pme_yz], rpole[k][mpl_pme_zz], f, 0, e, pgrad);


         if CONSTEXPR (do_a) {
            if (mscale == -1)
               atomic_add(-1, nem, offset);
         }
         if CONSTEXPR (do_e) {
            atomic_add(e, em, offset);
         }
         if CONSTEXPR (do_g) {
            atomic_add(pgrad.frcx, gx, i);
            atomic_add(pgrad.frcy, gy, i);
            atomic_add(pgrad.frcz, gz, i);
            atomic_add(-pgrad.frcx, gx, k);
            atomic_add(-pgrad.frcy, gy, k);
            atomic_add(-pgrad.frcz, gz, k);


            atomic_add(pgrad.ttmi[0], trqx, i);
            atomic_add(pgrad.ttmi[1], trqy, i);
            atomic_add(pgrad.ttmi[2], trqz, i);
            atomic_add(pgrad.ttmk[0], trqx, k);
            atomic_add(pgrad.ttmk[1], trqy, k);
            atomic_add(pgrad.ttmk[2], trqz, k);
         }
         if CONSTEXPR (do_v) {
            real vxx = -xr * pgrad.frcx;
            real vxy = -0.5f * (yr * pgrad.frcx + xr * pgrad.frcy);
            real vxz = -0.5f * (zr * pgrad.frcx + xr * pgrad.frcz);
            real vyy = -yr * pgrad.frcy;
            real vyz = -0.5f * (zr * pgrad.frcy + yr * pgrad.frcz);
            real vzz = -zr * pgrad.frcz;
            atomic_add(vxx, vxy, vxz, vyy, vyz, vzz, vir_em, offset);
         }
      } // end if (r2 <= off2)
   }
}


template <class Ver, class ETYP>
void empole_cu()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;


   const auto& st = *mspatial_unit;
   const real off = st.cutoff;
   const real off2 = off * off;
   auto bufsize = buffer_size();


   const real f = electric / dielec;
   real aewald = 0;
   if CONSTEXPR (eq<ETYP, EWALD>()) {
      PMEUnit pu = epme_unit;
      aewald = pu->aewald;


      if CONSTEXPR (do_e) {
         launch_k1s(nonblk, n, empole_self_cu<do_a>, //
                    bufsize, nem, em, rpole, n, f, aewald);
      }
   }
   if (st.niak > 0) {
      auto ker1 = empole_cu1<Ver, ETYP>;
      launch_k1s(nonblk, WARP_SIZE * st.niak, ker1, //
                 bufsize, nem, em, vir_em, gx, gy, gz, trqx, trqy, trqz,
                 TINKER_IMAGE_ARGS, off2, f, rpole, //
                 st.sorted, st.niak, st.iak, st.lst, n, aewald);
   }
   if (nmexclude_ > 0) {
      auto ker2 = empole_cu2<Ver>;
      launch_k1s(nonblk, nmexclude_, ker2, //
                 bufsize, nem, em, vir_em, gx, gy, gz, trqx, trqy, trqz,
                 TINKER_IMAGE_ARGS, off2, f, rpole, //
                 x, y, z, nmexclude_, mexclude_, mexclude_scale_);
   }
}


void empole_nonewald_cu(int vers)
{
   if (vers == calc::v0) {
      empole_cu<calc::V0, NON_EWALD>();
   } else if (vers == calc::v1) {
      empole_cu<calc::V1, NON_EWALD>();
   } else if (vers == calc::v3) {
      empole_cu<calc::V3, NON_EWALD>();
   } else if (vers == calc::v4) {
      empole_cu<calc::V4, NON_EWALD>();
   } else if (vers == calc::v5) {
      empole_cu<calc::V5, NON_EWALD>();
   } else if (vers == calc::v6) {
      empole_cu<calc::V6, NON_EWALD>();
   }
}


void empole_ewald_real_self_cu(int vers)
{
   if (vers == calc::v0) {
      empole_cu<calc::V0, EWALD>();
   } else if (vers == calc::v1) {
      empole_cu<calc::V1, EWALD>();
   } else if (vers == calc::v3) {
      empole_cu<calc::V3, EWALD>();
   } else if (vers == calc::v4) {
      empole_cu<calc::V4, EWALD>();
   } else if (vers == calc::v5) {
      empole_cu<calc::V5, EWALD>();
   } else if (vers == calc::v6) {
      empole_cu<calc::V6, EWALD>();
   }
}
TINKER_NAMESPACE_END
