#include "add.h"
#include "ff/amoeba/elecamoeba.h"
#include "ff/energy.h"
#include "ff/hippo/elechippo.h"
#include "ff/hippo/empolechgpen.h"
#include "ff/hippo/epolarchgpen.h"
#include "ff/image.h"
#include "ff/nblist.h"
#include "ff/pme.h"
#include "ff/potent.h"
#include "ff/switch.h"
#include "seq/pair_field_chgpen.h"
#include "tool/gpucard.h"

namespace tinker {
// see also subroutine udirect2b / dfield_chgpen0c in induce.f
#define DFIELD_DPTRS x, y, z, pcore, pval, palpha, field, rpole
template <class ETYP>
void dfield_chgpen_acc1(real (*field)[3])
{
   real aewald = 0;
   real off;
   if CONSTEXPR (eq<ETYP, EWALD>()) {
      off = switchOff(Switch::EWALD);
      const PMEUnit pu = ppme_unit;
      aewald = pu->aewald;
   } else {
      darray::zero(g::q0, n, field);
      off = switchOff(Switch::MPOLE);
   }

   const real off2 = off * off;
   const int maxnlst = mlist_unit->maxnlst;
   const auto* mlst = mlist_unit.deviceptr();

   MAYBE_UNUSED int GRID_DIM = gpuGridSize(BLOCK_DIM);
   #pragma acc parallel async num_gangs(GRID_DIM) vector_length(BLOCK_DIM)\
               present(lvec1,lvec2,lvec3,recipa,recipb,recipc)\
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
            pair_dfield_chgpen<ETYP>( //
               r2, xr, yr, zr, 1, ci, dix, diy, diz, corei, vali, alphai, qixx, qixy, qixz, qiyy,
               qiyz, qizz, rpole[k][mpl_pme_0], rpole[k][mpl_pme_x], rpole[k][mpl_pme_y],
               rpole[k][mpl_pme_z], pcore[k], pval[k], palpha[k], rpole[k][mpl_pme_xx],
               rpole[k][mpl_pme_xy], rpole[k][mpl_pme_xz], rpole[k][mpl_pme_yy],
               rpole[k][mpl_pme_yz], rpole[k][mpl_pme_zz], aewald, fid.x, fid.y, fid.z, fkd.x,
               fkd.y, fkd.z);

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
               present(lvec1,lvec2,lvec3,recipa,recipb,recipc)\
               deviceptr(DFIELD_DPTRS,mdwexclude,mdwexclude_scale)
   #pragma acc loop independent
   for (int ii = 0; ii < nmdwexclude; ++ii) {
      int i = mdwexclude[ii][0];
      int k = mdwexclude[ii][1];
      real dscale = mdwexclude_scale[ii][1] - 1;

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
         pair_dfield_chgpen<NON_EWALD>(r2, xr, yr, zr, dscale, ci, dix, diy, diz, corei, vali,
            alphai, qixx, qixy, qixz, qiyy, qiyz, qizz, rpole[k][mpl_pme_0], rpole[k][mpl_pme_x],
            rpole[k][mpl_pme_y], rpole[k][mpl_pme_z], pcore[k], pval[k], palpha[k],
            rpole[k][mpl_pme_xx], rpole[k][mpl_pme_xy], rpole[k][mpl_pme_xz], rpole[k][mpl_pme_yy],
            rpole[k][mpl_pme_yz], rpole[k][mpl_pme_zz], 0, fid.x, fid.y, fid.z, fkd.x, fkd.y,
            fkd.z);

         atomic_add(fid.x, &field[i][0]);
         atomic_add(fid.y, &field[i][1]);
         atomic_add(fid.z, &field[i][2]);

         atomic_add(fkd.x, &field[k][0]);
         atomic_add(fkd.y, &field[k][1]);
         atomic_add(fkd.z, &field[k][2]);
      }
   }
}

// see also subroutine umutual1 in induce.f
void ufield_chgpen_ewald_recip_self_acc(const real (*uind)[3], real (*field)[3])
{
   darray::zero(g::q0, n, field);

   const PMEUnit pu = ppme_unit;
   const auto& st = *pu;
   const int nfft1 = st.nfft1;
   const int nfft2 = st.nfft2;
   const int nfft3 = st.nfft3;
   const real aewald = st.aewald;

   cuind_to_fuind(pu, uind, uind, fuind, fuind);
   grid_uind(pu, fuind, fuind);
   fftfront(pu);
   // TODO: store vs. recompute qfac
   pme_conv(pu);
   fftback(pu);
   fphi_uind2(pu, fdip_phi1, fdip_phi2);

   const real term = aewald * aewald * aewald * 4 / 3 / sqrtpi;

   #pragma acc parallel loop independent async\
               present(lvec1,lvec2,lvec3,recipa,recipb,recipc)\
               deviceptr(field,uind,fdip_phi1)
   for (int i = 0; i < n; ++i) {
      real a[3][3];
      a[0][0] = nfft1 * recipa.x;
      a[1][0] = nfft2 * recipb.x;
      a[2][0] = nfft3 * recipc.x;
      a[0][1] = nfft1 * recipa.y;
      a[1][1] = nfft2 * recipb.y;
      a[2][1] = nfft3 * recipc.y;
      a[0][2] = nfft1 * recipa.z;
      a[1][2] = nfft2 * recipb.z;
      a[2][2] = nfft3 * recipc.z;

      #pragma acc loop seq
      for (int j = 0; j < 3; ++j) {
         real df1 =
            a[0][j] * fdip_phi1[i][1] + a[1][j] * fdip_phi1[i][2] + a[2][j] * fdip_phi1[i][3];
         field[i][j] += (term * uind[i][j] - df1);
      }
   }
}

#define UFIELD_DPTRS x, y, z, pcore, pval, palpha, field, uind
template <class ETYP>
void ufield_chgpen_acc1(const real (*uind)[3], real (*field)[3])
{
   real aewald = 0;
   real off;
   if CONSTEXPR (eq<ETYP, EWALD>()) {
      off = switchOff(Switch::EWALD);
      const PMEUnit pu = ppme_unit;
      aewald = pu->aewald;
   } else {
      darray::zero(g::q0, n, field);
      off = switchOff(Switch::MPOLE);
   }

   const real off2 = off * off;
   const int maxnlst = mlist_unit->maxnlst;
   const auto* mlst = mlist_unit.deviceptr();

   MAYBE_UNUSED int GRID_DIM = gpuGridSize(BLOCK_DIM);
   #pragma acc parallel async num_gangs(GRID_DIM) vector_length(BLOCK_DIM)\
               present(lvec1,lvec2,lvec3,recipa,recipb,recipc)\
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
            pair_ufield_chgpen<ETYP>(r2, xr, yr, zr, 1, uindi0, uindi1, uindi2, corei, vali, alphai,
               uind[k][0], uind[k][1], uind[k][2], pcore[k], pval[k], palpha[k], aewald, fid.x,
               fid.y, fid.z, fkd.x, fkd.y, fkd.z);

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
               present(lvec1,lvec2,lvec3,recipa,recipb,recipc)\
               deviceptr(UFIELD_DPTRS,mdwexclude,mdwexclude_scale)
   #pragma acc loop independent
   for (int ii = 0; ii < nmdwexclude; ++ii) {
      int i = mdwexclude[ii][0];
      int k = mdwexclude[ii][1];
      real wscale = mdwexclude_scale[ii][2] - 1;

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
         pair_ufield_chgpen<NON_EWALD>(r2, xr, yr, zr, wscale, uindi0, uindi1, uindi2, corei, vali,
            alphai, uind[k][0], uind[k][1], uind[k][2], pcore[k], pval[k], palpha[k], 0, fid.x,
            fid.y, fid.z, fkd.x, fkd.y, fkd.z);

         atomic_add(fid.x, &field[i][0]);
         atomic_add(fid.y, &field[i][1]);
         atomic_add(fid.z, &field[i][2]);

         atomic_add(fkd.x, &field[k][0]);
         atomic_add(fkd.y, &field[k][1]);
         atomic_add(fkd.z, &field[k][2]);
      }
   }
}

void dfield_chgpen_nonewald_acc(real (*field)[3])
{
   dfield_chgpen_acc1<NON_EWALD>(field);
}

void dfield_chgpen_ewald_real_acc(real (*field)[3])
{
   dfield_chgpen_acc1<EWALD>(field);
}

void ufield_chgpen_nonewald_acc(const real (*uind)[3], real (*field)[3])
{
   ufield_chgpen_acc1<NON_EWALD>(uind, field);
}

void ufield_chgpen_ewald_real_acc(const real (*uind)[3], real (*field)[3])
{
   ufield_chgpen_acc1<EWALD>(uind, field);
}
}
