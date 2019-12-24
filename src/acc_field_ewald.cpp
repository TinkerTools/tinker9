#include "acc_add.h"
#include "e_mpole.h"
#include "e_polar.h"
#include "gpu_card.h"
#include "md.h"
#include "nblist.h"
#include "pme.h"
#include "potent.h"
#include "seq_image.h"
#include "seq_pair_field.h"

TINKER_NAMESPACE_BEGIN
// see also subroutine udirect1 in induce.f
void dfield_ewald_recip_self(real (*field)[3])
{
   const PMEUnit pu = ppme_unit;
   const real aewald = pu->aewald;
   const real term = REAL_CUBE(aewald) * 4 / 3 / sqrtpi;

   cmp_to_fmp(pu, cmp, fmp);
   grid_mpole(pu, fmp);
   fftfront(pu);
   if (vir_m && !use_potent(mpole_term))
      pme_conv1(pu, vir_m);
   else
      pme_conv0(pu);
   fftback(pu);
   fphi_mpole(pu, fphi);
   fphi_to_cphi(pu, fphi, cphi);

   #pragma acc parallel loop independent async deviceptr(field,cphi,rpole)
   for (int i = 0; i < n; ++i) {
      real dix = rpole[i][mpl_pme_x];
      real diy = rpole[i][mpl_pme_y];
      real diz = rpole[i][mpl_pme_z];
      field[i][0] += (-cphi[i][1] + term * dix);
      field[i][1] += (-cphi[i][2] + term * diy);
      field[i][2] += (-cphi[i][3] + term * diz);
   }
}

// see also subroutine udirect2b / dfield0c in induce.f
void dfield_ewald_real_acc(real (*field)[3], real (*fieldp)[3])
{
   const real off = mlist_unit->cutoff;
   const real off2 = off * off;
   const int maxnlst = mlist_unit->maxnlst;
   const auto* mlst = mlist_unit.deviceptr();

   const PMEUnit pu = ppme_unit;
   const real aewald = pu->aewald;

#define DFIELD_DPTRS_ x, y, z, box, thole, pdamp, field, fieldp, rpole

   MAYBE_UNUSED int GRID_DIM = get_grid_size(BLOCK_DIM);
   #pragma acc parallel num_gangs(GRID_DIM) vector_length(BLOCK_DIM)\
               deviceptr(DFIELD_DPTRS_,mlst)
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
      real pdi = pdamp[i];
      real pti = thole[i];
      real gxi = 0, gyi = 0, gzi = 0;
      real txi = 0, tyi = 0, tzi = 0;

      int nmlsti = mlst->nlst[i];
      int base = i * maxnlst;
      #pragma acc loop vector independent reduction(+:gxi,gyi,gzi,txi,tyi,tzi)
      for (int kk = 0; kk < nmlsti; ++kk) {
         int k = mlst->lst[base + kk];
         real xr = x[k] - xi;
         real yr = y[k] - yi;
         real zr = z[k] - zi;

         image(xr, yr, zr, box);
         real r2 = xr * xr + yr * yr + zr * zr;
         if (r2 <= off2) {
            real3 fid = make_real3(0, 0, 0);
            real3 fip = make_real3(0, 0, 0);
            real3 fkd = make_real3(0, 0, 0);
            real3 fkp = make_real3(0, 0, 0);
            pair_dfield<elec_t::ewald>( //
               r2, xr, yr, zr, 1, 1,    //
               ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz, pdi,
               pti, //
               rpole[k][mpl_pme_0], rpole[k][mpl_pme_x], rpole[k][mpl_pme_y],
               rpole[k][mpl_pme_z], rpole[k][mpl_pme_xx], rpole[k][mpl_pme_xy],
               rpole[k][mpl_pme_xz], rpole[k][mpl_pme_yy], rpole[k][mpl_pme_yz],
               rpole[k][mpl_pme_zz], pdamp[k], thole[k], //
               aewald, fid, fip, fkd, fkp);

            gxi += fid.x;
            gyi += fid.y;
            gzi += fid.z;
            txi += fip.x;
            tyi += fip.y;
            tzi += fip.z;

            atomic_add(fkd.x, &field[k][0]);
            atomic_add(fkd.y, &field[k][1]);
            atomic_add(fkd.z, &field[k][2]);
            atomic_add(fkp.x, &fieldp[k][0]);
            atomic_add(fkp.y, &fieldp[k][1]);
            atomic_add(fkp.z, &fieldp[k][2]);
         }
      } // end for (int kk)

      atomic_add(gxi, &field[i][0]);
      atomic_add(gyi, &field[i][1]);
      atomic_add(gzi, &field[i][2]);
      atomic_add(txi, &fieldp[i][0]);
      atomic_add(tyi, &fieldp[i][1]);
      atomic_add(tzi, &fieldp[i][2]);
   } // end for (int i)

   #pragma acc parallel deviceptr(DFIELD_DPTRS_,dpexclude_,dpexclude_scale_)
   #pragma acc loop independent
   for (int ii = 0; ii < ndpexclude_; ++ii) {
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
         real3 fid = make_real3(0, 0, 0);
         real3 fip = make_real3(0, 0, 0);
         real3 fkd = make_real3(0, 0, 0);
         real3 fkp = make_real3(0, 0, 0);
         pair_dfield<elec_t::coulomb>(                                       //
            r2, xr, yr, zr, dscale, pscale,                                  //
            ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz, pdi, pti, //
            rpole[k][mpl_pme_0], rpole[k][mpl_pme_x], rpole[k][mpl_pme_y],
            rpole[k][mpl_pme_z], rpole[k][mpl_pme_xx], rpole[k][mpl_pme_xy],
            rpole[k][mpl_pme_xz], rpole[k][mpl_pme_yy], rpole[k][mpl_pme_yz],
            rpole[k][mpl_pme_zz], pdamp[k], thole[k], //
            0, fid, fip, fkd, fkp);

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
      }
   }
}

// see also subroutine umutual1 in induce.f
void ufield_ewald_recip_self(const real (*uind)[3], const real (*uinp)[3],
                             real (*field)[3], real (*fieldp)[3])
{
   const PMEUnit pu = ppme_unit;
   const auto& st = *pu;
   const int nfft1 = st.nfft1;
   const int nfft2 = st.nfft2;
   const int nfft3 = st.nfft3;
   const real aewald = st.aewald;

   cuind_to_fuind(pu, uind, uinp, fuind, fuinp);
   grid_uind(pu, fuind, fuinp);
   fftfront(pu);
   // TODO: store vs. recompute qfac
   pme_conv0(pu);
   fftback(pu);
   fphi_uind2(pu, fdip_phi1, fdip_phi2);

   const real term = REAL_CUBE(aewald) * 4 / 3 / sqrtpi;

   #pragma acc parallel loop independent async\
               deviceptr(box,field,fieldp,uind,uinp,fdip_phi1,fdip_phi2)
   for (int i = 0; i < n; ++i) {
      real a[3][3];
      a[0][0] = nfft1 * box->recip[0][0];
      a[1][0] = nfft2 * box->recip[1][0];
      a[2][0] = nfft3 * box->recip[2][0];
      a[0][1] = nfft1 * box->recip[0][1];
      a[1][1] = nfft2 * box->recip[1][1];
      a[2][1] = nfft3 * box->recip[2][1];
      a[0][2] = nfft1 * box->recip[0][2];
      a[1][2] = nfft2 * box->recip[1][2];
      a[2][2] = nfft3 * box->recip[2][2];

      #pragma acc loop seq
      for (int j = 0; j < 3; ++j) {
         real df1 = a[0][j] * fdip_phi1[i][1] + a[1][j] * fdip_phi1[i][2] +
            a[2][j] * fdip_phi1[i][3];
         real df2 = a[0][j] * fdip_phi2[i][1] + a[1][j] * fdip_phi2[i][2] +
            a[2][j] * fdip_phi2[i][3];
         field[i][j] += (term * uind[i][j] - df1);
         fieldp[i][j] += (term * uinp[i][j] - df2);
      }
   }
}

void ufield_ewald_real_acc(const real (*uind)[3], const real (*uinp)[3],
                           real (*field)[3], real (*fieldp)[3])
{
   const real off = mlist_unit->cutoff;
   const real off2 = off * off;
   const int maxnlst = mlist_unit->maxnlst;
   const auto* mlst = mlist_unit.deviceptr();

   const PMEUnit pu = ppme_unit;
   const real aewald = pu->aewald;

#define UFIELD_DPTRS_ x, y, z, box, thole, pdamp, field, fieldp, uind, uinp

   MAYBE_UNUSED int GRID_DIM = get_grid_size(BLOCK_DIM);
   #pragma acc parallel num_gangs(GRID_DIM) vector_length(BLOCK_DIM)\
               deviceptr(UFIELD_DPTRS_,mlst)
   #pragma acc loop gang independent
   for (int i = 0; i < n; ++i) {
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real uindi0 = uind[i][0];
      real uindi1 = uind[i][1];
      real uindi2 = uind[i][2];
      real uinpi0 = uinp[i][0];
      real uinpi1 = uinp[i][1];
      real uinpi2 = uinp[i][2];
      real pdi = pdamp[i];
      real pti = thole[i];
      real gxi = 0, gyi = 0, gzi = 0;
      real txi = 0, tyi = 0, tzi = 0;

      int nmlsti = mlst->nlst[i];
      int base = i * maxnlst;
      #pragma acc loop vector independent reduction(+:gxi,gyi,gzi,txi,tyi,tzi)
      for (int kk = 0; kk < nmlsti; ++kk) {
         int k = mlst->lst[base + kk];
         real xr = x[k] - xi;
         real yr = y[k] - yi;
         real zr = z[k] - zi;

         image(xr, yr, zr, box);
         real r2 = xr * xr + yr * yr + zr * zr;
         if (r2 <= off2) {
            real3 fid = make_real3(0, 0, 0);
            real3 fip = make_real3(0, 0, 0);
            real3 fkd = make_real3(0, 0, 0);
            real3 fkp = make_real3(0, 0, 0);
            pair_ufield<elec_t::ewald>(                                  //
               r2, xr, yr, zr, 1,                                        //
               uindi0, uindi1, uindi2, uinpi0, uinpi1, uinpi2, pdi, pti, //
               uind[k][0], uind[k][1], uind[k][2], uinp[k][0], uinp[k][1],
               uinp[k][2], pdamp[k], thole[k], //
               aewald, fid, fip, fkd, fkp);

            gxi += fid.x;
            gyi += fid.y;
            gzi += fid.z;
            txi += fip.x;
            tyi += fip.y;
            tzi += fip.z;

            atomic_add(fkd.x, &field[k][0]);
            atomic_add(fkd.y, &field[k][1]);
            atomic_add(fkd.z, &field[k][2]);
            atomic_add(fkp.x, &fieldp[k][0]);
            atomic_add(fkp.y, &fieldp[k][1]);
            atomic_add(fkp.z, &fieldp[k][2]);
         }
      } // end for (int kk)

      atomic_add(gxi, &field[i][0]);
      atomic_add(gyi, &field[i][1]);
      atomic_add(gzi, &field[i][2]);
      atomic_add(txi, &fieldp[i][0]);
      atomic_add(tyi, &fieldp[i][1]);
      atomic_add(tzi, &fieldp[i][2]);
   } // end for (int i)

   #pragma acc parallel deviceptr(UFIELD_DPTRS_,uexclude_,uexclude_scale_)
   #pragma acc loop independent
   for (int ii = 0; ii < nuexclude_; ++ii) {
      int i = uexclude_[ii][0];
      int k = uexclude_[ii][1];
      real uscale = uexclude_scale_[ii];

      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real uindi0 = uind[i][0];
      real uindi1 = uind[i][1];
      real uindi2 = uind[i][2];
      real uinpi0 = uinp[i][0];
      real uinpi1 = uinp[i][1];
      real uinpi2 = uinp[i][2];
      real pdi = pdamp[i];
      real pti = thole[i];

      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;

      image(xr, yr, zr, box);
      real r2 = xr * xr + yr * yr + zr * zr;
      if (r2 <= off2) {
         real3 fid = make_real3(0, 0, 0);
         real3 fip = make_real3(0, 0, 0);
         real3 fkd = make_real3(0, 0, 0);
         real3 fkp = make_real3(0, 0, 0);
         pair_ufield<elec_t::coulomb>(
            r2, xr, yr, zr, uscale,                                   //
            uindi0, uindi1, uindi2, uinpi0, uinpi1, uinpi2, pdi, pti, //
            uind[k][0], uind[k][1], uind[k][2], uinp[k][0], uinp[k][1],
            uinp[k][2], pdamp[k], thole[k], //
            0, fid, fip, fkd, fkp);

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
      }
   }
}
TINKER_NAMESPACE_END
