#include "add.cuh"
#include "e_polar.h"
#include "launch.cuh"
#include "md.h"
#include "pme.h"
#include "seq_image.h"
#include "seq_pair_field.h"
#include "spatial.h"
#include "switch.h"


TINKER_NAMESPACE_BEGIN
#define DFIELD_ARGS                                                            \
   real(*restrict field)[3], real(*restrict fieldp)[3],                        \
      const real *restrict thole, const real *restrict pdamp,                  \
      const real(*restrict rpole)[10], real aewald, const Box *restrict box,   \
      real off
#define SHFL_SYMB(s) __shfl_sync(ALL_LANES, sh##s, srclane)


__global__
void dfield_ewald_real_cu1(DFIELD_ARGS, const Spatial* restrict sp)
{
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);


   real gxi, gyi, gzi, txi, tyi, tzi;
   __shared__ real gxk[BLOCK_DIM], gyk[BLOCK_DIM], gzk[BLOCK_DIM],
      txk[BLOCK_DIM], tyk[BLOCK_DIM], tzk[BLOCK_DIM];


   const real off2 = off * off;
   const int n = sp->n;
   const int niak = sp->niak;
   const auto* restrict sorted = sp->sorted;
   const auto* restrict iak = sp->iak;
   const auto* restrict lst = sp->lst;
   for (int iw = iwarp; iw < niak; iw += nwarp) {
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
      real shpdk = pdamp[shk];
      real shptk = thole[shk];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int atomk = SHFL_SYMB(atomk);
         real xr = __shfl_sync(ALL_LANES, shx, srclane) - xi;
         real yr = __shfl_sync(ALL_LANES, shy, srclane) - yi;
         real zr = __shfl_sync(ALL_LANES, shz, srclane) - zi;
         int k = SHFL_SYMB(k);
         real ck = SHFL_SYMB(ck);
         real dkx = SHFL_SYMB(dkx);
         real dky = SHFL_SYMB(dky);
         real dkz = SHFL_SYMB(dkz);
         real qkxx = SHFL_SYMB(qkxx);
         real qkxy = SHFL_SYMB(qkxy);
         real qkxz = SHFL_SYMB(qkxz);
         real qkyy = SHFL_SYMB(qkyy);
         real qkyz = SHFL_SYMB(qkyz);
         real qkzz = SHFL_SYMB(qkzz);
         real pdk = SHFL_SYMB(pdk);
         real ptk = SHFL_SYMB(ptk);


         PairField pairf;
         pairf.fid[0] = 0;
         pairf.fid[1] = 0;
         pairf.fid[2] = 0;
         pairf.fkd[0] = 0;
         pairf.fkd[1] = 0;
         pairf.fkd[2] = 0;
         pairf.fip[0] = 0;
         pairf.fip[1] = 0;
         pairf.fip[2] = 0;
         pairf.fkp[0] = 0;
         pairf.fkp[1] = 0;
         pairf.fkp[2] = 0;


         image(xr, yr, zr, box);
         real r2 = xr * xr + yr * yr + zr * zr;
         if (atomi < atomk && r2 <= off2) {
            pair_dfield<elec_t::ewald>(
               r2, xr, yr, zr, 1, 1, ci, dix, diy, diz, qixx, qixy, qixz, qiyy,
               qiyz, qizz, pdi, pti, ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy,
               qkyz, qkzz, pdk, ptk, aewald, pairf);
         } // end if (r2 <= off2)


         gxi += pairf.fid[0];
         gyi += pairf.fid[1];
         gzi += pairf.fid[2];
         gxk[srclane + (threadIdx.x - ilane)] += pairf.fkd[0];
         gyk[srclane + (threadIdx.x - ilane)] += pairf.fkd[1];
         gzk[srclane + (threadIdx.x - ilane)] += pairf.fkd[2];
         txi += pairf.fip[0];
         tyi += pairf.fip[1];
         tzi += pairf.fip[2];
         txk[srclane + (threadIdx.x - ilane)] += pairf.fkp[0];
         tyk[srclane + (threadIdx.x - ilane)] += pairf.fkp[1];
         tzk[srclane + (threadIdx.x - ilane)] += pairf.fkp[2];
      } // end for (j)


      atomic_add(gxi, &field[i][0]);
      atomic_add(gyi, &field[i][1]);
      atomic_add(gzi, &field[i][2]);
      atomic_add(txi, &fieldp[i][0]);
      atomic_add(tyi, &fieldp[i][1]);
      atomic_add(tzi, &fieldp[i][2]);
      atomic_add(gxk[threadIdx.x], &field[shk][0]);
      atomic_add(gyk[threadIdx.x], &field[shk][1]);
      atomic_add(gzk[threadIdx.x], &field[shk][2]);
      atomic_add(txk[threadIdx.x], &fieldp[shk][0]);
      atomic_add(tyk[threadIdx.x], &fieldp[shk][1]);
      atomic_add(tzk[threadIdx.x], &fieldp[shk][2]);
   } // end for (iw)
}


__global__
void dfield_real_cu2(DFIELD_ARGS, const real* restrict x,
                     const real* restrict y, const real* restrict z,
                     int ndpexclude_, const int (*restrict dpexclude_)[2],
                     const real (*restrict dpexclude_scale_)[2])
{
   const real off2 = off * off;
   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < ndpexclude_;
        ii += blockDim.x * gridDim.x) {
      int i = dpexclude_[ii][0];
      int k = dpexclude_[ii][1];
      real dscale = dpexclude_scale_[ii][0];
      real pscale = dpexclude_scale_[ii][1];


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
      real pdi = pdamp[i];
      real pti = thole[i];


      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;


      image(xr, yr, zr, box);
      real r2 = xr * xr + yr * yr + zr * zr;
      if (r2 <= off2) {
         PairField pairf;
         pair_dfield<elec_t::coulomb>(
            r2, xr, yr, zr, dscale, pscale, ci, dix, diy, diz, qixx, qixy, qixz,
            qiyy, qiyz, qizz, pdi, pti, rpole[k][mpl_pme_0],
            rpole[k][mpl_pme_x], rpole[k][mpl_pme_y], rpole[k][mpl_pme_z],
            rpole[k][mpl_pme_xx], rpole[k][mpl_pme_xy], rpole[k][mpl_pme_xz],
            rpole[k][mpl_pme_yy], rpole[k][mpl_pme_yz], rpole[k][mpl_pme_zz],
            pdamp[k], thole[k], 0, pairf);


         atomic_add(pairf.fid[0], &field[i][0]);
         atomic_add(pairf.fid[1], &field[i][1]);
         atomic_add(pairf.fid[2], &field[i][2]);
         atomic_add(pairf.fip[0], &fieldp[i][0]);
         atomic_add(pairf.fip[1], &fieldp[i][1]);
         atomic_add(pairf.fip[2], &fieldp[i][2]);


         atomic_add(pairf.fkd[0], &field[k][0]);
         atomic_add(pairf.fkd[1], &field[k][1]);
         atomic_add(pairf.fkd[2], &field[k][2]);
         atomic_add(pairf.fkp[0], &fieldp[k][0]);
         atomic_add(pairf.fkp[1], &fieldp[k][1]);
         atomic_add(pairf.fkp[2], &fieldp[k][2]);
      } // end if (r2 <= off2)
   }    // end for (ii)
}


void dfield_ewald_real_cu(real (*field)[3], real (*fieldp)[3])
{
   const real off = switch_off(switch_ewald);
   const auto& st = *mspatial_unit;
   const auto* sp = mspatial_unit.deviceptr();


   const PMEUnit pu = ppme_unit;
   const real aewald = pu->aewald;


   launch_kernel1(WARP_SIZE * st.niak, dfield_ewald_real_cu1, field, fieldp,
                  thole, pdamp, rpole, aewald, box, off, sp);
   if (ndpexclude_ > 0)
      launch_kernel1(ndpexclude_, dfield_real_cu2, field, fieldp, thole, pdamp,
                     rpole, aewald, box, off, x, y, z, ndpexclude_, dpexclude_,
                     dpexclude_scale_);
}
TINKER_NAMESPACE_END
