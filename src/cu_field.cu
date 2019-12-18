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
      const real(*restrict rpole)[10], const Box *restrict box, real off


template <elec_t ETYP>
__global__
void dfield_cu1(DFIELD_ARGS, const Spatial* restrict sp, real aewald)
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
         real pdk = __shfl_sync(ALL_LANES, shpdk, srclane);
         real ptk = __shfl_sync(ALL_LANES, shptk, srclane);


         PairField pairf;
         zero(pairf);


         image(xr, yr, zr, box);
         real r2 = xr * xr + yr * yr + zr * zr;
         if (atomi < atomk && r2 <= off2) {
            if CONSTEXPR (ETYP == elec_t::ewald) {
               pair_dfield<elec_t::ewald>(
                  r2, xr, yr, zr, 1, 1, ci, dix, diy, diz, qixx, qixy, qixz,
                  qiyy, qiyz, qizz, pdi, pti, ck, dkx, dky, dkz, qkxx, qkxy,
                  qkxz, qkyy, qkyz, qkzz, pdk, ptk, aewald, pairf);
            }
            if CONSTEXPR (ETYP == elec_t::coulomb) {
               pair_dfield<elec_t::coulomb>(
                  r2, xr, yr, zr, 1, 1, ci, dix, diy, diz, qixx, qixy, qixz,
                  qiyy, qiyz, qizz, pdi, pti, ck, dkx, dky, dkz, qkxx, qkxy,
                  qkxz, qkyy, qkyz, qkzz, pdk, ptk, 0, pairf);
            }
         } // end if (include)


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
void dfield_cu2(DFIELD_ARGS, const real* restrict x, const real* restrict y,
                const real* restrict z, int ndpexclude_,
                const int (*restrict dpexclude_)[2],
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
      } // end if (include)
   }    // end for (ii)
}


void dfield_ewald_real_cu(real (*field)[3], real (*fieldp)[3])
{
   const auto& st = *mspatial_unit;
   const real off = st.cutoff;
   const auto* sp = mspatial_unit.deviceptr();


   const PMEUnit pu = ppme_unit;
   const real aewald = pu->aewald;


   if (st.niak > 0) {
      launch_kernel1(WARP_SIZE * st.niak, dfield_cu1<elec_t::ewald>, field,
                     fieldp, thole, pdamp, rpole, box, off, sp, aewald);
   }
   if (ndpexclude_ > 0) {
      launch_kernel1(ndpexclude_, dfield_cu2, field, fieldp, thole, pdamp,
                     rpole, box, off, x, y, z, ndpexclude_, dpexclude_,
                     dpexclude_scale_);
   }
}


void dfield_coulomb_cu(real (*field)[3], real (*fieldp)[3])
{
   const auto& st = *mspatial_unit;
   const real off = st.cutoff;
   const auto* sp = mspatial_unit.deviceptr();


   device_array::zero(n, field, fieldp);
   if (st.niak > 0) {
      launch_kernel1(WARP_SIZE * st.niak, dfield_cu1<elec_t::coulomb>, field,
                     fieldp, thole, pdamp, rpole, box, off, sp, 0);
   }
   if (ndpexclude_ > 0) {
      launch_kernel1(ndpexclude_, dfield_cu2, field, fieldp, thole, pdamp,
                     rpole, box, off, x, y, z, ndpexclude_, dpexclude_,
                     dpexclude_scale_);
   }
}


#define UFIELD_ARGS                                                            \
   const real(*restrict uind)[3], const real(*restrict uinp)[3],               \
      real(*restrict field)[3], real(*restrict fieldp)[3],                     \
      const real *restrict thole, const real *restrict pdamp,                  \
      TINKER_IMAGE_PARAMS, real off2


template <elec_t ETYP>
__launch_bounds__(BLOCK_DIM) __global__
void ufield_cu1(UFIELD_ARGS, int n, const Spatial::SortedAtom* restrict sorted,
                int niak, const int* restrict iak, const int* restrict lst,
                real aewald)
{
   const int iwarp = (threadIdx.x + blockIdx.x * blockDim.x) / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);


   struct Data
   {
      real3 fkd, fkp, ukd, ukp, rk;
      real pdk, ptk;
   };
   __shared__ Data data[BLOCK_DIM];


   for (int iw = iwarp; iw < niak; iw += nwarp) {
      real3 fid = make_real3(0, 0, 0);
      real3 fip = make_real3(0, 0, 0);
      int atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
      real3 ri = make_real3(sorted[atomi].x, sorted[atomi].y, sorted[atomi].z);
      int i = sorted[atomi].unsorted;
      real3 uid = make_real3(uind[i][0], uind[i][1], uind[i][2]);
      real3 uip = make_real3(uinp[i][0], uinp[i][1], uinp[i][2]);
      real pdi = pdamp[i];
      real pti = thole[i];


      data[threadIdx.x].fkd = make_real3(0, 0, 0);
      data[threadIdx.x].fkp = make_real3(0, 0, 0);
      int shatomk = lst[iw * WARP_SIZE + ilane];
      data[threadIdx.x].rk =
         make_real3(sorted[shatomk].x, sorted[shatomk].y, sorted[shatomk].z);
      int shk = sorted[shatomk].unsorted;
      data[threadIdx.x].ukd =
         make_real3(uind[shk][0], uind[shk][1], uind[shk][2]);
      data[threadIdx.x].ukp =
         make_real3(uinp[shk][0], uinp[shk][1], uinp[shk][2]);
      data[threadIdx.x].pdk = pdamp[shk];
      data[threadIdx.x].ptk = thole[shk];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
         real3 dr = data[klane].rk - ri;


         real r2 = image2(dr.x, dr.y, dr.z);
         if (atomi < atomk && r2 <= off2) {
            if CONSTEXPR (ETYP == elec_t::ewald) {
               pair_ufield<elec_t::ewald>(
                  r2, dr.x, dr.y, dr.z, 1, uid.x, uid.y, uid.z, uip.x, uip.y,
                  uip.z, pdi, pti, data[klane].ukd.x, data[klane].ukd.y,
                  data[klane].ukd.z, data[klane].ukp.x, data[klane].ukp.y,
                  data[klane].ukp.z, data[klane].pdk, data[klane].ptk, aewald,
                  fid, fip, data[klane].fkd, data[klane].fkp);
            }
            if CONSTEXPR (ETYP == elec_t::coulomb) {
               pair_ufield<elec_t::coulomb>(
                  r2, dr.x, dr.y, dr.z, 1, uid.x, uid.y, uid.z, uip.x, uip.y,
                  uip.z, pdi, pti, data[klane].ukd.x, data[klane].ukd.y,
                  data[klane].ukd.z, data[klane].ukp.x, data[klane].ukp.y,
                  data[klane].ukp.z, data[klane].pdk, data[klane].ptk, 0, fid,
                  fip, data[klane].fkd, data[klane].fkp);
            }
         } // end if (include)
      }    // end for (j)


      atomic_add(fid.x, &field[i][0]);
      atomic_add(fid.y, &field[i][1]);
      atomic_add(fid.z, &field[i][2]);
      atomic_add(fip.x, &fieldp[i][0]);
      atomic_add(fip.y, &fieldp[i][1]);
      atomic_add(fip.z, &fieldp[i][2]);
      atomic_add(data[threadIdx.x].fkd.x, &field[shk][0]);
      atomic_add(data[threadIdx.x].fkd.y, &field[shk][1]);
      atomic_add(data[threadIdx.x].fkd.z, &field[shk][2]);
      atomic_add(data[threadIdx.x].fkp.x, &fieldp[shk][0]);
      atomic_add(data[threadIdx.x].fkp.y, &fieldp[shk][1]);
      atomic_add(data[threadIdx.x].fkp.z, &fieldp[shk][2]);
   } // end for (iw)
}


__global__
void ufield_cu2(UFIELD_ARGS, const real* restrict x, const real* restrict y,
                const real* restrict z, int nuexclude_,
                const int (*restrict uexclude_)[2],
                const real* restrict uexclude_scale_)
{
   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < nuexclude_;
        ii += blockDim.x * gridDim.x) {
      int i = uexclude_[ii][0];
      int k = uexclude_[ii][1];
      real uscale = uexclude_scale_[ii];


      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real3 uid = make_real3(uind[i][0], uind[i][1], uind[i][2]);
      real3 uip = make_real3(uinp[i][0], uinp[i][1], uinp[i][2]);
      real pdi = pdamp[i];
      real pti = thole[i];


      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;


      real r2 = image2(xr, yr, zr);
      if (r2 <= off2) {
         real3 fid = make_real3(0, 0, 0);
         real3 fip = make_real3(0, 0, 0);
         real3 fkd = make_real3(0, 0, 0);
         real3 fkp = make_real3(0, 0, 0);
         pair_ufield<elec_t::coulomb>(
            r2, xr, yr, zr, uscale, uid.x, uid.y, uid.z, uip.x, uip.y, uip.z,
            pdi, pti, uind[k][0], uind[k][1], uind[k][2], uinp[k][0],
            uinp[k][1], uinp[k][2], pdamp[k], thole[k], 0, fid, fip, fkd, fkp);


         atomic_add(fid.x, &field[i][0]);
         atomic_add(fid.y, &field[i][1]);
         atomic_add(fid.z, &field[i][2]);
         atomic_add(fip.x, &fieldp[i][0]);
         atomic_add(fip.y, &fieldp[i][1]);
         atomic_add(fip.z, &fieldp[i][2]);


         atomic_add(fkd.x, &field[k][0]);
         atomic_add(fkd.y, &field[k][1]);
         atomic_add(fkd.z, &field[k][2]);
         atomic_add(fkp.x, &fieldp[k][0]);
         atomic_add(fkp.y, &fieldp[k][1]);
         atomic_add(fkp.z, &fieldp[k][2]);
      } // end if (include)
   }
}


void ufield_ewald_real_cu(const real (*uind)[3], const real (*uinp)[3],
                          real (*field)[3], real (*fieldp)[3])
{
   const auto& st = *mspatial_unit;
   const real off2 = st.cutoff * st.cutoff;


   const PMEUnit pu = ppme_unit;
   const real aewald = pu->aewald;

   if (st.niak > 0) {
      launch_kernel1(WARP_SIZE * st.niak, ufield_cu1<elec_t::ewald>, uind, uinp,
                     field, fieldp, thole, pdamp, TINKER_IMAGE_ARGS, off2, n,
                     st.sorted, st.niak, st.iak, st.lst, aewald);
   }
   if (nuexclude_) {
      launch_kernel1(nuexclude_, ufield_cu2, uind, uinp, field, fieldp, thole,
                     pdamp, TINKER_IMAGE_ARGS, off2, x, y, z, nuexclude_,
                     uexclude_, uexclude_scale_);
   }
}


void ufield_coulomb_cu(float const (*uind)[3], float const (*uinp)[3],
                       float (*field)[3], float (*fieldp)[3])
{
   const auto& st = *mspatial_unit;
   const real off2 = st.cutoff * st.cutoff;


   device_array::zero(n, field, fieldp);
   if (st.niak > 0) {
      launch_kernel1(WARP_SIZE * st.niak, ufield_cu1<elec_t::coulomb>, uind,
                     uinp, field, fieldp, thole, pdamp, TINKER_IMAGE_ARGS, off2,
                     n, st.sorted, st.niak, st.iak, st.lst, 0);
   }
   if (nuexclude_) {
      launch_kernel1(nuexclude_, ufield_cu2, uind, uinp, field, fieldp, thole,
                     pdamp, TINKER_IMAGE_ARGS, off2, x, y, z, nuexclude_,
                     uexclude_, uexclude_scale_);
   }
}
TINKER_NAMESPACE_END
