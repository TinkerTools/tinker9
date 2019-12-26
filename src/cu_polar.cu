#include "add.h"
#include "e_polar.h"
#include "epolar_trq.h"
#include "launch.h"
#include "md.h"
#include "pme.h"
#include "seq_image.h"
#include "seq_pair_polar.h"
#include "spatial.h"


TINKER_NAMESPACE_BEGIN
#define POLAR_ARGS                                                             \
   size_t bufsize, count_buffer restrict nep, energy_buffer restrict ep,       \
      virial_buffer restrict vir_ep, real *restrict gx, real *restrict gy,     \
      real *restrict gz, real(*restrict ufld)[3], real(*restrict dufld)[6],    \
      const Box *restrict box, real off2, real f,                              \
      const real(*restrict rpole)[10], const real *restrict pdamp,             \
      const real *restrict thole, const real(*restrict uind)[3],               \
      const real(*restrict uinp)[3]


template <int USE, elec_t ETYP>
__global__
void epolar_cu1(POLAR_ARGS, const Spatial::SortedAtom* restrict sorted,
                int niak, const int* restrict iak, const int* restrict lst,
                int n, real aewald)
{
   constexpr int do_e = USE & calc::energy;
   constexpr int do_a = USE & calc::analyz;
   constexpr int do_g = USE & calc::grad;
   constexpr int do_v = USE & calc::virial;


   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);
   const int offset = ithread & (bufsize - 1);


   MAYBE_UNUSED int ctl;
   MAYBE_UNUSED real etl;
   MAYBE_UNUSED real vtlxx, vtlxy, vtlxz, vtlyy, vtlyz, vtlzz;
   MAYBE_UNUSED real gxi, gyi, gzi, txi, tyi, tzi, dui[6];
   MAYBE_UNUSED __shared__ real gxk[BLOCK_DIM], gyk[BLOCK_DIM], gzk[BLOCK_DIM],
      txk[BLOCK_DIM], tyk[BLOCK_DIM], tzk[BLOCK_DIM], duk[BLOCK_DIM][6];


   for (int iw = iwarp; iw < niak; iw += nwarp) {
      if CONSTEXPR (do_a && do_e) {
         ctl = 0;
         etl = 0;
      }
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
         #pragma unroll
         for (int i = 0; i < 6; ++i) {
            dui[i] = 0;
         }
         gxk[threadIdx.x] = 0;
         gyk[threadIdx.x] = 0;
         gzk[threadIdx.x] = 0;
         txk[threadIdx.x] = 0;
         tyk[threadIdx.x] = 0;
         tzk[threadIdx.x] = 0;
         #pragma unroll
         for (int i = 0; i < 6; ++i) {
            duk[threadIdx.x][i] = 0;
         }
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
      real uix = uind[i][0];
      real uiy = uind[i][1];
      real uiz = uind[i][2];
      real uixp = uinp[i][0];
      real uiyp = uinp[i][1];
      real uizp = uinp[i][2];
      real pdi = pdamp[i];
      real pti = thole[i];


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
      real shukx = uind[shk][0];
      real shuky = uind[shk][1];
      real shukz = uind[shk][2];
      real shukxp = uinp[shk][0];
      real shukyp = uinp[shk][1];
      real shukzp = uinp[shk][2];
      real shpdk = pdamp[shk];
      real shptk = thole[shk];


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
         real ukx = __shfl_sync(ALL_LANES, shukx, srclane);
         real uky = __shfl_sync(ALL_LANES, shuky, srclane);
         real ukz = __shfl_sync(ALL_LANES, shukz, srclane);
         real ukxp = __shfl_sync(ALL_LANES, shukxp, srclane);
         real ukyp = __shfl_sync(ALL_LANES, shukyp, srclane);
         real ukzp = __shfl_sync(ALL_LANES, shukzp, srclane);
         real pdk = __shfl_sync(ALL_LANES, shpdk, srclane);
         real ptk = __shfl_sync(ALL_LANES, shptk, srclane);


         real e;
         PairPolarGrad pgrad;
         zero(pgrad);


         image(xr, yr, zr, box);
         real r2 = xr * xr + yr * yr + zr * zr;
         if (atomi < atomk && r2 <= off2) {
            if CONSTEXPR (ETYP == elec_t::ewald) {
               pair_polar<USE, elec_t::ewald>( //
                  r2, xr, yr, zr, 1, 1, 1,     //
                  ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz, uix,
                  uiy, uiz, uixp, uiyp, uizp, pdi, pti, //
                  ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukx,
                  uky, ukz, ukxp, ukyp, ukzp, pdk, ptk, //
                  f, aewald, e, pgrad);
            }
            if CONSTEXPR (ETYP == elec_t::coulomb) {
               pair_polar<USE, elec_t::coulomb>( //
                  r2, xr, yr, zr, 1, 1, 1,       //
                  ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz, uix,
                  uiy, uiz, uixp, uiyp, uizp, pdi, pti, //
                  ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukx,
                  uky, ukz, ukxp, ukyp, ukzp, pdk, ptk, //
                  f, 0, e, pgrad);
            }


            if CONSTEXPR (do_a && do_e) {
               ctl += 1;
               etl += e;
            }
            if CONSTEXPR (do_v) {
               vtlxx += -xr * pgrad.frcx;
               vtlxy += -0.5f * (yr * pgrad.frcx + xr * pgrad.frcy);
               vtlxz += -0.5f * (zr * pgrad.frcx + xr * pgrad.frcz);
               vtlyy += -yr * pgrad.frcy;
               vtlyz += -0.5f * (zr * pgrad.frcy + yr * pgrad.frcz);
               vtlzz += -zr * pgrad.frcz;
            }
         } // end if (include)


         if CONSTEXPR (do_g) {
            gxi += pgrad.frcx;
            gyi += pgrad.frcy;
            gzi += pgrad.frcz;
            gxk[srclane + (threadIdx.x - ilane)] -= pgrad.frcx;
            gyk[srclane + (threadIdx.x - ilane)] -= pgrad.frcy;
            gzk[srclane + (threadIdx.x - ilane)] -= pgrad.frcz;


            txi += pgrad.ufldi[0];
            tyi += pgrad.ufldi[1];
            tzi += pgrad.ufldi[2];
            txk[srclane + (threadIdx.x - ilane)] += pgrad.ufldk[0];
            tyk[srclane + (threadIdx.x - ilane)] += pgrad.ufldk[1];
            tzk[srclane + (threadIdx.x - ilane)] += pgrad.ufldk[2];


            dui[0] += pgrad.dufldi[0];
            dui[1] += pgrad.dufldi[1];
            dui[2] += pgrad.dufldi[2];
            dui[3] += pgrad.dufldi[3];
            dui[4] += pgrad.dufldi[4];
            dui[5] += pgrad.dufldi[5];
            duk[srclane + (threadIdx.x - ilane)][0] += pgrad.dufldk[0];
            duk[srclane + (threadIdx.x - ilane)][1] += pgrad.dufldk[1];
            duk[srclane + (threadIdx.x - ilane)][2] += pgrad.dufldk[2];
            duk[srclane + (threadIdx.x - ilane)][3] += pgrad.dufldk[3];
            duk[srclane + (threadIdx.x - ilane)][4] += pgrad.dufldk[4];
            duk[srclane + (threadIdx.x - ilane)][5] += pgrad.dufldk[5];
         }
      } // end for (j)


      if CONSTEXPR (do_a && do_e) {
         atomic_add(ctl, nep, offset);
         atomic_add(etl, ep, offset);
      }
      if CONSTEXPR (do_g) {
         atomic_add(gxi, gx, i);
         atomic_add(gyi, gy, i);
         atomic_add(gzi, gz, i);
         atomic_add(txi, &ufld[i][0]);
         atomic_add(tyi, &ufld[i][1]);
         atomic_add(tzi, &ufld[i][2]);
         atomic_add(gxk[threadIdx.x], gx, shk);
         atomic_add(gyk[threadIdx.x], gy, shk);
         atomic_add(gzk[threadIdx.x], gz, shk);
         atomic_add(txk[threadIdx.x], &ufld[shk][0]);
         atomic_add(tyk[threadIdx.x], &ufld[shk][1]);
         atomic_add(tzk[threadIdx.x], &ufld[shk][2]);


         atomic_add(dui[0], &dufld[i][0]);
         atomic_add(dui[1], &dufld[i][1]);
         atomic_add(dui[2], &dufld[i][2]);
         atomic_add(dui[3], &dufld[i][3]);
         atomic_add(dui[4], &dufld[i][4]);
         atomic_add(dui[5], &dufld[i][5]);
         atomic_add(duk[threadIdx.x][0], &dufld[shk][0]);
         atomic_add(duk[threadIdx.x][1], &dufld[shk][1]);
         atomic_add(duk[threadIdx.x][2], &dufld[shk][2]);
         atomic_add(duk[threadIdx.x][3], &dufld[shk][3]);
         atomic_add(duk[threadIdx.x][4], &dufld[shk][4]);
         atomic_add(duk[threadIdx.x][5], &dufld[shk][5]);
      }
      if CONSTEXPR (do_v)
         atomic_add(vtlxx, vtlxy, vtlxz, vtlyy, vtlyz, vtlzz, vir_ep, offset);
   } // end for (iw)
}


template <int USE>
__global__
void epolar_cu2(POLAR_ARGS, const real* restrict x, const real* restrict y,
                const real* restrict z, int ndpuexclude_,
                const int (*restrict dpuexclude_)[2],
                const real (*restrict dpuexclude_scale_)[3])
{
   constexpr int do_e = USE & calc::energy;
   constexpr int do_a = USE & calc::analyz;
   constexpr int do_g = USE & calc::grad;
   constexpr int do_v = USE & calc::virial;


   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < ndpuexclude_;
        ii += blockDim.x * gridDim.x) {
      int offset = ii & (bufsize - 1);


      int i = dpuexclude_[ii][0];
      int k = dpuexclude_[ii][1];
      real dscale = dpuexclude_scale_[ii][0];
      real pscale = dpuexclude_scale_[ii][1];
      real uscale = dpuexclude_scale_[ii][2];


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
      real uix = uind[i][0];
      real uiy = uind[i][1];
      real uiz = uind[i][2];
      real pdi = pdamp[i];
      real pti = thole[i];
      MAYBE_UNUSED real uixp, uiyp, uizp;
      if CONSTEXPR (do_g) {
         uixp = uinp[i][0];
         uiyp = uinp[i][1];
         uizp = uinp[i][2];
      }


      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;


      image(xr, yr, zr, box);
      real r2 = xr * xr + yr * yr + zr * zr;
      if (r2 <= off2) {
         real ck = rpole[k][mpl_pme_0];
         real dkx = rpole[k][mpl_pme_x];
         real dky = rpole[k][mpl_pme_y];
         real dkz = rpole[k][mpl_pme_z];
         real qkxx = rpole[k][mpl_pme_xx];
         real qkxy = rpole[k][mpl_pme_xy];
         real qkxz = rpole[k][mpl_pme_xz];
         real qkyy = rpole[k][mpl_pme_yy];
         real qkyz = rpole[k][mpl_pme_yz];
         real qkzz = rpole[k][mpl_pme_zz];
         real ukx = uind[k][0];
         real uky = uind[k][1];
         real ukz = uind[k][2];
         MAYBE_UNUSED real ukxp, ukyp, ukzp;
         if CONSTEXPR (do_g) {
            ukxp = uinp[k][0];
            ukyp = uinp[k][1];
            ukzp = uinp[k][2];
         }


         real e;
         PairPolarGrad pgrad;
         pair_polar<USE, elec_t::coulomb>(          //
            r2, xr, yr, zr, dscale, pscale, uscale, //
            ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz, uix, uiy,
            uiz, uixp, uiyp, uizp, pdi, pti, //
            ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukx, uky,
            ukz, ukxp, ukyp, ukzp, pdamp[k], thole[k], //
            f, 0, e, pgrad);


         if CONSTEXPR (do_a && do_e) {
            if (pscale == -1)
               atomic_add(-1, nep, offset);
            atomic_add(e, ep, offset);
         }
         if CONSTEXPR (do_g) {
            atomic_add(pgrad.frcx, gx, i);
            atomic_add(pgrad.frcy, gy, i);
            atomic_add(pgrad.frcz, gz, i);
            atomic_add(-pgrad.frcx, gx, k);
            atomic_add(-pgrad.frcy, gy, k);
            atomic_add(-pgrad.frcz, gz, k);


            atomic_add(pgrad.ufldi[0], &ufld[i][0]);
            atomic_add(pgrad.ufldi[1], &ufld[i][1]);
            atomic_add(pgrad.ufldi[2], &ufld[i][2]);
            atomic_add(pgrad.ufldk[0], &ufld[k][0]);
            atomic_add(pgrad.ufldk[1], &ufld[k][1]);
            atomic_add(pgrad.ufldk[2], &ufld[k][2]);


            atomic_add(pgrad.dufldi[0], &dufld[i][0]);
            atomic_add(pgrad.dufldi[1], &dufld[i][1]);
            atomic_add(pgrad.dufldi[2], &dufld[i][2]);
            atomic_add(pgrad.dufldi[3], &dufld[i][3]);
            atomic_add(pgrad.dufldi[4], &dufld[i][4]);
            atomic_add(pgrad.dufldi[5], &dufld[i][5]);
            atomic_add(pgrad.dufldk[0], &dufld[k][0]);
            atomic_add(pgrad.dufldk[1], &dufld[k][1]);
            atomic_add(pgrad.dufldk[2], &dufld[k][2]);
            atomic_add(pgrad.dufldk[3], &dufld[k][3]);
            atomic_add(pgrad.dufldk[4], &dufld[k][4]);
            atomic_add(pgrad.dufldk[5], &dufld[k][5]);


            if CONSTEXPR (do_v) {
               real vxx = -xr * pgrad.frcx;
               real vxy = -0.5f * (yr * pgrad.frcx + xr * pgrad.frcy);
               real vxz = -0.5f * (zr * pgrad.frcx + xr * pgrad.frcz);
               real vyy = -yr * pgrad.frcy;
               real vyz = -0.5f * (zr * pgrad.frcy + yr * pgrad.frcz);
               real vzz = -zr * pgrad.frcz;
               atomic_add(vxx, vxy, vxz, vyy, vyz, vzz, vir_ep, offset);
            }
         }
      } // end if (r2 <= off2)
   }
}


template <int USE, elec_t ETYP>
void epolar_tmpl_cu(const real (*uind)[3], const real (*uinp)[3])
{
   constexpr int do_g = USE & calc::grad;


   const auto& st = *mspatial_unit;
   const real off = st.cutoff;
   const real off2 = off * off;
   auto bufsize = buffer_size();


   const real f = 0.5f * electric / dielec;
   real aewald = 0;
   if CONSTEXPR (ETYP == elec_t::ewald) {
      PMEUnit pu = ppme_unit;
      aewald = pu->aewald;
   }


   if CONSTEXPR (do_g) {
      device_array::zero(n, ufld, dufld);
   }
   if (st.niak > 0) {
      auto ker1 = epolar_cu1<USE, ETYP>;
      launch_k1(WARP_SIZE * st.niak, ker1, //
                bufsize, nep, ep, vir_ep, gx, gy, gz, ufld, dufld, box, off2, f,
                rpole, pdamp, thole, uind, uinp, //
                st.sorted, st.niak, st.iak, st.lst, n, aewald);
   }
   if (ndpuexclude_ > 0) {
      auto ker2 = epolar_cu2<USE>;
      launch_k1(ndpuexclude_, ker2, //
                bufsize, nep, ep, vir_ep, gx, gy, gz, ufld, dufld, box, off2, f,
                rpole, pdamp, thole, uind, uinp, //
                x, y, z, ndpuexclude_, dpuexclude_, dpuexclude_scale_);
   }
   // torque
   if CONSTEXPR (do_g) {
      launch_k1(n, epolar_trq_cu, //
                trqx, trqy, trqz, n, rpole, ufld, dufld);
   }
}


template <int USE>
void epolar_real_cu(const real (*uind)[3], const real (*uinp)[3])
{
   epolar_tmpl_cu<USE, elec_t::ewald>(uind, uinp);
}
template <>
void epolar_real_cu<calc::v0>(const real (*uind)[3], const real (*)[3])
{
   epolar0_dotprod(uind, udirp);
}
template void epolar_real_cu<calc::v1>(const real (*)[3], const real (*)[3]);
template void epolar_real_cu<calc::v3>(const real (*)[3], const real (*)[3]);
template void epolar_real_cu<calc::v4>(const real (*)[3], const real (*)[3]);
template void epolar_real_cu<calc::v5>(const real (*)[3], const real (*)[3]);
template void epolar_real_cu<calc::v6>(const real (*)[3], const real (*)[3]);


void epolar_coulomb_cu(int vers, const real (*uind)[3], const real (*uinp)[3])
{
   if (vers == calc::v0) {
      epolar0_dotprod(uind, udirp);
   } else if (vers == calc::v1) {
      epolar_tmpl_cu<calc::v1, elec_t::coulomb>(uind, uinp);
   } else if (vers == calc::v3) {
      epolar_tmpl_cu<calc::v3, elec_t::coulomb>(uind, uinp);
   } else if (vers == calc::v4) {
      epolar_tmpl_cu<calc::v4, elec_t::coulomb>(uind, uinp);
   } else if (vers == calc::v5) {
      epolar_tmpl_cu<calc::v5, elec_t::coulomb>(uind, uinp);
   } else if (vers == calc::v6) {
      epolar_tmpl_cu<calc::v6, elec_t::coulomb>(uind, uinp);
   }
}
TINKER_NAMESPACE_END
