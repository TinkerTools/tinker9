#include "add.h"
#include "empole_chgpen.h"
#include "epolar_chgpen.h"
#include "glob.nblist.h"
#include "image.h"
#include "md.h"
#include "seq_pair_field_chgpen.h"
#include "switch.h"
#include "tool/gpu_card.h"

namespace tinker {
// see also subroutine dfield_chgpen0b in induce.f
#define DFIELD_DPTRS x, y, z, pcore, pval, palpha, field, rpole
// TODO: HIPPO not reviewed
void dfield_chgpen_nonewald_acc(real (*field)[3])
{
   darray::zero(PROCEED_NEW_Q, n, field);

   const real off = switch_off(switch_mpole);
   const real off2 = off * off;
   const int maxnlst = mlist_unit->maxnlst;
   const auto* mlst = mlist_unit.deviceptr();

   MAYBE_UNUSED int GRID_DIM = get_grid_size(BLOCK_DIM);
   #pragma acc parallel async num_gangs(GRID_DIM) vector_length(BLOCK_DIM)\
               deviceptr(DFIELD_DPTRS,mlst)
   #pragma acc loop gang independent
   for (int i = 0; i < n; ++i) {
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
      real gxi = 0, gyi = 0, gzi = 0;

      int nmlsti = mlst->nlst[i];
      int base = i * maxnlst;
      #pragma acc loop vector independent reduction(+:gxi,gyi,gzi)
      for (int kk = 0; kk < nmlsti; ++kk) {
         int k = mlst->lst[base + kk];
         real xr = x[k] - xi;
         real yr = y[k] - yi;
         real zr = z[k] - zi;

         real r2 = image2(xr, yr, zr);
         if (r2 <= off2) {
            real3 fid = make_real3(0, 0, 0);
            real3 fkd = make_real3(0, 0, 0);
            pair_dfield_chgpen<NON_EWALD>( //
               r2, xr, yr, zr, 1, ci, dix, diy, diz, corei, vali, alphai, qixx,
               qixy, qixz, qiyy, qiyz, qizz, rpole[k][mpl_pme_0],
               rpole[k][mpl_pme_x], rpole[k][mpl_pme_y], rpole[k][mpl_pme_z],
               pcore[k], pval[k], palpha[k], rpole[k][mpl_pme_xx],
               rpole[k][mpl_pme_xy], rpole[k][mpl_pme_xz], rpole[k][mpl_pme_yy],
               rpole[k][mpl_pme_yz], rpole[k][mpl_pme_zz], 0, fid.x, fid.y,
               fid.z, fkd.x, fkd.y, fkd.z);

            gxi += fid.x;
            gyi += fid.y;
            gzi += fid.z;


            atomic_add(fkd.x, &field[k][0]);
            atomic_add(fkd.y, &field[k][1]);
            atomic_add(fkd.z, &field[k][2]);
         }
      } // end for (int kk)

      atomic_add(gxi, &field[i][0]);
      atomic_add(gyi, &field[i][1]);
      atomic_add(gzi, &field[i][2]);
   } // end for (int i)

   #pragma acc parallel async\
               deviceptr(DFIELD_DPTRS,dexclude,dexclude_scale)
   #pragma acc loop independent
   for (int ii = 0; ii < ndexclude; ++ii) {
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
         pair_dfield_chgpen<NON_EWALD>( //
            r2, xr, yr, zr, dscale, ci, dix, diy, diz, corei, vali, alphai,
            qixx, qixy, qixz, qiyy, qiyz, qizz, rpole[k][mpl_pme_0],
            rpole[k][mpl_pme_x], rpole[k][mpl_pme_y], rpole[k][mpl_pme_z],
            pcore[k], pval[k], palpha[k], rpole[k][mpl_pme_xx],
            rpole[k][mpl_pme_xy], rpole[k][mpl_pme_xz], rpole[k][mpl_pme_yy],
            rpole[k][mpl_pme_yz], rpole[k][mpl_pme_zz], 0, fid.x, fid.y, fid.z,
            fkd.x, fkd.y, fkd.z);

         atomic_add(fid.x, &field[i][0]);
         atomic_add(fid.y, &field[i][1]);
         atomic_add(fid.z, &field[i][2]);

         atomic_add(fkd.x, &field[k][0]);
         atomic_add(fkd.y, &field[k][1]);
         atomic_add(fkd.z, &field[k][2]);
      }
   }
}


// see also subroutine ufield0b in induce.f
#define UFIELD_DPTRS x, y, z, pcore, pval, palpha, field, uind
// TODO: HIPPO not reviewed
void ufield_chgpen_nonewald_acc(const real (*uind)[3], real (*field)[3])
{
   darray::zero(PROCEED_NEW_Q, n, field);

   const real off = switch_off(switch_mpole);
   const real off2 = off * off;
   const int maxnlst = mlist_unit->maxnlst;
   const auto* mlst = mlist_unit.deviceptr();

   MAYBE_UNUSED int GRID_DIM = get_grid_size(BLOCK_DIM);
   #pragma acc parallel async num_gangs(GRID_DIM) vector_length(BLOCK_DIM)\
               deviceptr(UFIELD_DPTRS,mlst)
   #pragma acc loop gang independent
   for (int i = 0; i < n; ++i) {
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real uindi0 = uind[i][0];
      real uindi1 = uind[i][1];
      real uindi2 = uind[i][2];
      real corei = pcore[i];
      real alphai = palpha[i];
      real vali = pval[i];
      real gxi = 0, gyi = 0, gzi = 0;

      int nmlsti = mlst->nlst[i];
      int base = i * maxnlst;
      #pragma acc loop vector independent reduction(+:gxi,gyi,gzi)
      for (int kk = 0; kk < nmlsti; ++kk) {
         int k = mlst->lst[base + kk];
         real xr = x[k] - xi;
         real yr = y[k] - yi;
         real zr = z[k] - zi;

         real r2 = image2(xr, yr, zr);
         if (r2 <= off2) {
            real3 fid = make_real3(0, 0, 0);
            real3 fkd = make_real3(0, 0, 0);
            pair_ufield_chgpen<NON_EWALD>(
               r2, xr, yr, zr, 1, uindi0, uindi1, uindi2, corei, vali, alphai,
               uind[k][0], uind[k][1], uind[k][2], pcore[k], pval[k], palpha[k],
               0, fid.x, fid.y, fid.z, fkd.x, fkd.y, fkd.z);

            gxi += fid.x;
            gyi += fid.y;
            gzi += fid.z;


            atomic_add(fkd.x, &field[k][0]);
            atomic_add(fkd.y, &field[k][1]);
            atomic_add(fkd.z, &field[k][2]);
         }
      } // end for (int kk)

      atomic_add(gxi, &field[i][0]);
      atomic_add(gyi, &field[i][1]);
      atomic_add(gzi, &field[i][2]);

   } // end for (int i)

   #pragma acc parallel async\
               deviceptr(UFIELD_DPTRS,wexclude,wexclude_scale)
   #pragma acc loop independent
   for (int ii = 0; ii < nwexclude; ++ii) {
      int i = wexclude[ii][0];
      int k = wexclude[ii][1];
      real wscale = wexclude_scale[ii] - 1;

      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real uindi0 = uind[i][0];
      real uindi1 = uind[i][1];
      real uindi2 = uind[i][2];
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
            r2, xr, yr, zr, wscale, uindi0, uindi1, uindi2, corei, vali, alphai,
            uind[k][0], uind[k][1], uind[k][2], pcore[k], pval[k], palpha[k], 0,
            fid.x, fid.y, fid.z, fkd.x, fkd.y, fkd.z);

         atomic_add(fid.x, &field[i][0]);
         atomic_add(fid.y, &field[i][1]);
         atomic_add(fid.z, &field[i][2]);

         atomic_add(fkd.x, &field[k][0]);
         atomic_add(fkd.y, &field[k][1]);
         atomic_add(fkd.z, &field[k][2]);
      }
   }
}
}
