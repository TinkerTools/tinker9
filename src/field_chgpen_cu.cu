#include "add.h"
#include "epolar_chgpen.h"
#include "glob.spatial.h"
#include "image.h"
#include "launch.h"
#include "md.h"
#include "pme.h"
#include "seq_pair_field_chgpen.h"
#include "switch.h"
#include "tool/cudalib.h"


namespace tinker {
#define DFIELDPARAS                                                            \
   real(*restrict field)[3], const real *restrict pcore,                       \
      const real *restrict pval, const real *restrict palpha,
const real (*restrict rpole)[10], TINKER_IMAGE_PARAMS,
   real off2


   template <class ETYP>
   __launch_bounds__(BLOCK_DIM) __global__
void dfield_chgpen_cu1(DFIELDPARAS, int n,
                       const Spatial::SortedAtom* restrict sorted, int niak,
                       const int* restrict iak, const int* restrict lst,
                       real aewald)
{
   const int iwarp = (threadIdx.x + blockIdx.x * blockDim.x) / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);


   struct Data
   {
      real3 fkd;
      real3 rk;
      real ck;
      real dkx, dky, dkz;
      real qkxx, qkxy, qkxz, qkyy, qkyz, qkzz;
      real core, val, alpha;
   };
   __shared__ Data data[BLOCK_DIM];


   for (int iw = iwarp; iw < niak; iw += nwarp) {
      real3 fid = make_real3(0, 0, 0);
      int atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
      real3 ri = make_real3(sorted[atomi].x, sorted[atomi].y, sorted[atomi].z);
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
      real corei = pcore[i];
      real alphai = palpha[i];
      real vali = pval[i];


      data[threadIdx.x].fkd = make_real3(0, 0, 0);
      data[threadIdx.x].fkp = make_real3(0, 0, 0);
      int shatomk = lst[iw * WARP_SIZE + ilane];
      data[threadIdx.x].rk =
         make_real3(sorted[shatomk].x, sorted[shatomk].y, sorted[shatomk].z);
      int shk = sorted[shatomk].unsorted;
      data[threadIdx.x].ck = rpole[shk][mpl_pme_0];
      data[threadIdx.x].dkx = rpole[shk][mpl_pme_x];
      data[threadIdx.x].dky = rpole[shk][mpl_pme_y];
      data[threadIdx.x].dkz = rpole[shk][mpl_pme_z];
      data[threadIdx.x].qkxx = rpole[shk][mpl_pme_xx];
      data[threadIdx.x].qkxy = rpole[shk][mpl_pme_xy];
      data[threadIdx.x].qkxz = rpole[shk][mpl_pme_xz];
      data[threadIdx.x].qkyy = rpole[shk][mpl_pme_yy];
      data[threadIdx.x].qkyz = rpole[shk][mpl_pme_yz];
      data[threadIdx.x].qkzz = rpole[shk][mpl_pme_zz];
      data[threadIdx.x].core = pcore[shk];
      data[threadIdx.x].val = pval[shk];
      data[threadIdx.x].alpha = palpha[shk];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
         real3 dr = data[klane].rk - ri;


         real r2 = image2(dr.x, dr.y, dr.z);
         if (atomi < atomk && r2 <= off2) {
            if CONSTEXPR (eq<ETYP, EWALD>()) {
               pair_dfield_chgpen<EWALD>(
                  r2, dr.x, dr.y, dr.z, 1, ci, dix, diy, diz, corei, vali,
                  alphai, qixx, qixy, qixz, qiyy, qiyz, qizz, data[klane].ck,
                  data[klane].dkx, data[klane].dky, data[klane].dkz,
                  data[klane].core, data[klane].val, data[klane].alpha,
                  data[klane].qkxx, data[klane].qkxy, data[klane].qkxz,
                  data[klane].qkyy, data[klane].qkyz, data[klane].qkzz,
                  data[klane].pdk, data[klane].ptk, aewald, fid,
                  data[klane].fkd);
            }
            if CONSTEXPR (eq<ETYP, NON_EWALD>()) {
               pair_dfield_chgpen<NON_EWALD>(
                  r2, dr.x, dr.y, dr.z, 1, ci, dix, diy, diz, corei, vali,
                  alphai, qixx, qixy, qixz, qiyy, qiyz, qizz, data[klane].ck,
                  data[klane].dkx, data[klane].dky, data[klane].dkz,
                  data[klane].core, data[klane].val, data[klane].alpha,
                  data[klane].qkxx, data[klane].qkxy, data[klane].qkxz,
                  data[klane].qkyy, data[klane].qkyz, data[klane].qkzz,
                  data[klane].pdk, data[klane].ptk, 0, fid,
                  data[klane].fkd);
            }
         } // end if (include)
      }


      atomic_add(fid.x, &field[i][0]);
      atomic_add(fid.y, &field[i][1]);
      atomic_add(fid.z, &field[i][2]);
      atomic_add(data[threadIdx.x].fkd.x, &field[shk][0]);
      atomic_add(data[threadIdx.x].fkd.y, &field[shk][1]);
      atomic_add(data[threadIdx.x].fkd.z, &field[shk][2]);
   } // end for (iw)
}


__global__
void dfield_chgpen_cu2(DFIELDPARAS, const real* restrict x,
                       const real* restrict y, const real* restrict z,
                       int ndexclude, const int (*restrict dexclude)[2],
                       const real (*restrict dexclude_scale)[2])
{
   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < ndexclude;
        ii += blockDim.x * gridDim.x) {
      int i = dexclude[ii][0];
      int k = dexclude[ii][1];
      real dscale = dexclude_scale[ii];


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
      real corei = pcore[i];
      real alphai = palpha[i];
      real vali = pval[i];


      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;


      real r2 = image2(xr, yr, zr);
      if (r2 <= off2) {
         real3 fid = make_real3(0, 0, 0);
         real3 fkd = make_real3(0, 0, 0);
         pair_dfield_chgpen<NON_EWALD>(
            r2, dr.x, dr.y, dr.z, dscale, ci, dix, diy, diz, corei, vali,
            alphai, qixx, qixy, qixz, qiyy, qiyz, qizz, rpole[k][mpl_pme_0],
            rpole[k][mpl_pme_x], rpole[k][mpl_pme_y], rpole[k][mpl_pme_z],
            pcore[k], pval[k], palpha[k], rpole[k][mpl_pme_xx],
            rpole[k][mpl_pme_xy], rpole[k][mpl_pme_xz], rpole[k][mpl_pme_yy],
            rpole[k][mpl_pme_yz], rpole[k][mpl_pme_zz], 0, fid, fkd);


         atomic_add(fid.x, &field[i][0]);
         atomic_add(fid.y, &field[i][1]);
         atomic_add(fid.z, &field[i][2]);


         atomic_add(fkd.x, &field[k][0]);
         atomic_add(fkd.y, &field[k][1]);
         atomic_add(fkd.z, &field[k][2]);
      } // end if (include)
   }
}


void dfield_chgpen_ewald_real_cu(real (*field)[3])
{
   const auto& st = *mspatial_unit;
   const real off = switch_off(switch_ewald);
   const real off2 = off * off;


   const PMEUnit pu = ppme_unit;
   const real aewald = pu->aewald;
   if (st.niak > 0) {
      launch_k1s(nonblk, WARP_SIZE * st.niak, dfield_chgpen_cu1<EWALD>, field,
                 pcore, pval, palpha, rpole, TINKER_IMAGE_ARGS, off2, n,
                 st.sorted, st.niak, st.iak, st.lst, aewald);
   }
   if (ndexclude > 0) {
      launch_k1s(nonblk, ndexclude, dfield_chgpen_cu2, field, pcore, pval,
                 palpha, rpole, TINKER_IMAGE_ARGS, off2, x, y, z, ndexclude,
                 dexclude, dexclude_scale);
   }
}


void dfield_chgpen_nonewald_cu(real (*field)[3])
{
   const auto& st = *mspatial_unit;
   const real off = switch_off(switch_mpole);
   const real off2 = off * off;


   darray::zero(PROCEED_NEW_Q, n, field);
   if (st.niak > 0) {
      launch_k1s(nonblk, WARP_SIZE * st.niak, dfield_chgpen_cu1<NON_EWALD>,
                 field, pcore, pval, palpha, rpole, TINKER_IMAGE_ARGS, off2, n,
                 st.sorted, st.niak, st.iak, st.lst, 0);
   }
   if (ndexclude > 0) {
      launch_k1s(nonblk, ndexclude, dfield_chgpen_cu2, field, pcore, pval,
                 palpha, rpole, TINKER_IMAGE_ARGS, off2, x, y, z, ndexclude,
                 dexclude, dexclude_scale);
   }
}


#define UFIELDPARAS                                                            \
   const real(*restrict uind)[3], real(*restrict field)[3],                    \
      const real *restrict pcore, const real *restrict pval,                   \
      const real *restrict palpha,
TINKER_IMAGE_PARAMS,
   real off2


   template <class ETYP>
   __launch_bounds__(BLOCK_DIM) __global__
void ufield_chgpen_cu1(UFIELDPARAS, int n,
                       const Spatial::SortedAtom* restrict sorted, int niak,
                       const int* restrict iak, const int* restrict lst,
                       real aewald)
{
   const int iwarp = (threadIdx.x + blockIdx.x * blockDim.x) / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);


   struct Data
   {
      real3 fkd, ukd, rk;
      real core, alpha, val;
   };
   __shared__ Data data[BLOCK_DIM];


   for (int iw = iwarp; iw < niak; iw += nwarp) {
      real3 fid = make_real3(0, 0, 0);
      int atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
      real3 ri = make_real3(sorted[atomi].x, sorted[atomi].y, sorted[atomi].z);
      int i = sorted[atomi].unsorted;
      real3 uid = make_real3(uind[i][0], uind[i][1], uind[i][2]);
      real corei = pcore[i];
      real alphai = palpha[i];
      real vali = pval[i];


      data[threadIdx.x].fkd = make_real3(0, 0, 0);
      int shatomk = lst[iw * WARP_SIZE + ilane];
      data[threadIdx.x].rk =
         make_real3(sorted[shatomk].x, sorted[shatomk].y, sorted[shatomk].z);
      int shk = sorted[shatomk].unsorted;
      data[threadIdx.x].ukd =
         make_real3(uind[shk][0], uind[shk][1], uind[shk][2]);
      data[threadIdx.x].core = pcore[shk];
      data[threadIdx.x].val = pval[shk];
      data[threadIdx.x].alpha = palpha[shk];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
         real3 dr = data[klane].rk - ri;


         real r2 = image2(dr.x, dr.y, dr.z);
         if (atomi < atomk && r2 <= off2) {
            if CONSTEXPR (eq<ETYP, EWALD>()) {
               pair_ufield_chgpen<EWALD>(
                  r2, dr.x, dr.y, dr.z, 1, uid.x, uid.y, uid.z, corei, vali,
                  alphai, data[klane].ukd.x, data[klane].ukd.y,
                  data[klane].ukd.z, data[klane].core, data[klane].val,
                  data[klane].alpha, aewald, fid, data[klane].fkd);
            }
            if CONSTEXPR (eq<ETYP, NON_EWALD>()) {
               pair_ufield_chgpen<NON_EWALD>(
                  r2, dr.x, dr.y, dr.z, 1, uid.x, uid.y, uid.z, corei, vali,
                  alphai, data[klane].ukd.x, data[klane].ukd.y,
                  data[klane].ukd.z, data[klane].core, data[klane].val,
                  data[klane].alpha, aewald, fid, data[klane].fkd);
            }
         } // end if (include)
      }


      atomic_add(fid.x, &field[i][0]);
      atomic_add(fid.y, &field[i][1]);
      atomic_add(fid.z, &field[i][2]);
      atomic_add(data[threadIdx.x].fkd.x, &field[shk][0]);
      atomic_add(data[threadIdx.x].fkd.y, &field[shk][1]);
      atomic_add(data[threadIdx.x].fkd.z, &field[shk][2]);
   } // end for (iw)
}


__global__
void ufield_chgpen_cu2(UFIELDPARAS, const real* restrict x,
                       const real* restrict y, const real* restrict z,
                       int nwexclude, const int (*restrict wexclude)[2],
                       const real* restrict wexclude_scale)
{
   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < nwexclude;
        ii += blockDim.x * gridDim.x) {
      int i = wexclude[ii][0];
      int k = wexclude[ii][1];
      real wscale = wexclude_scale[ii];


      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real3 uid = make_real3(uind[i][0], uind[i][1], uind[i][2]);
      real corei = pcore[i];
      real alphai = palpha[i];
      real vali = pval[i];


      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;


      real r2 = image2(xr, yr, zr);
      if (r2 <= off2) {
         real3 fid = make_real3(0, 0, 0);
         real3 fkd = make_real3(0, 0, 0);
         pair_ufield_chgpen<NON_EWALD>(
            r2, xr, yr, zr, wscale, uid.x, uid.y, uid.z, corei, vali,
            alphai, uind[k][0], uind[k][1], uind[k][2],
            pcore[k], pval[k], palpha[k], 0, fid, fkd);


         atomic_add(fid.x, &field[i][0]);
         atomic_add(fid.y, &field[i][1]);
         atomic_add(fid.z, &field[i][2]);


         atomic_add(fkd.x, &field[k][0]);
         atomic_add(fkd.y, &field[k][1]);
         atomic_add(fkd.z, &field[k][2]);
      } // end if (include)
   }
}


void ufield_chgpen_ewald_real_cu(const real (*uind)[3], real (*field)[3])
{
   const auto& st = *mspatial_unit;
   const real off = switch_off(switch_ewald);
   const real off2 = off * off;


   const PMEUnit pu = ppme_unit;
   const real aewald = pu->aewald;


   if (st.niak > 0) {
      launch_k1s(nonblk, WARP_SIZE * st.niak, ufield_chgpen_cu1<EWALD>, uind,
                 field, pcore, pval, palpha, TINKER_IMAGE_ARGS, off2, n,
                 st.sorted, st.niak, st.iak, st.lst, aewald);
   }
   if (nwexclude) {
      launch_k1s(nonblk, nwexclude, ufield_chgpen_cu2, uind, field, pcore, pval,
                 palpha, TINKER_IMAGE_ARGS, off2, x, y, z, nwexclude, wexclude,
                 wexclude_scale);
   }
}


void ufield_chgpen_nonewald_cu(const real (*uind)[3], real (*field)[3])
{
   const auto& st = *mspatial_unit;
   const real off = switch_off(switch_mpole);
   const real off2 = off * off;


   darray::zero(PROCEED_NEW_Q, n, field);
   if (st.niak > 0) {
      launch_k1s(nonblk, WARP_SIZE * st.niak, ufield_chgpen_cu1<NON_EWALD>,
                 uind, field, pcore, pval, palpha, TINKER_IMAGE_ARGS, off2, n,
                 st.sorted, st.niak, st.iak, st.lst, 0);
   }
   if (nwexclude) {
      launch_k1s(nonblk, nwexclude, ufield_chgpen_cu2, uind, field, pcore, pval,
                 palpha, TINKER_IMAGE_ARGS, off2, x, y, z, nwexclude, wexclude,
                 wexclude_scale);
   }
}
}
