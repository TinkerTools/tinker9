#include "acc_add.h"
#include "e_polar.h"
#include "gpu_card.h"
#include "md.h"
#include "nblist.h"
#include "seq_image.h"
#include "seq_pair_field.h"
#include "seq_switch.h"

TINKER_NAMESPACE_BEGIN
// see also subroutine dfield0b in induce.f
void dfield_coulomb(real (*field)[3], real (*fieldp)[3])
{
   device_array::zero(n, field, fieldp);

   const real off = switch_off(switch_mpole);
   const real off2 = off * off;
   const int maxnlst = mlist_unit->maxnlst;
   const auto* mlst = mlist_unit.deviceptr();

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
            PairField pairf;
            pair_dfield<elec_t::coulomb>( //
               r2, xr, yr, zr, 1, 1,      //
               ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz, pdi,
               pti, //
               rpole[k][mpl_pme_0], rpole[k][mpl_pme_x], rpole[k][mpl_pme_y],
               rpole[k][mpl_pme_z], rpole[k][mpl_pme_xx], rpole[k][mpl_pme_xy],
               rpole[k][mpl_pme_xz], rpole[k][mpl_pme_yy], rpole[k][mpl_pme_yz],
               rpole[k][mpl_pme_zz], pdamp[k], thole[k], //
               0, pairf);

            gxi += pairf.fid[0];
            gyi += pairf.fid[1];
            gzi += pairf.fid[2];
            txi += pairf.fip[0];
            tyi += pairf.fip[1];
            tzi += pairf.fip[2];

            atomic_add(pairf.fkd[0], &field[k][0]);
            atomic_add(pairf.fkd[1], &field[k][1]);
            atomic_add(pairf.fkd[2], &field[k][2]);
            atomic_add(pairf.fkp[0], &fieldp[k][0]);
            atomic_add(pairf.fkp[1], &fieldp[k][1]);
            atomic_add(pairf.fkp[2], &fieldp[k][2]);
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

      PairField pairf;
      pair_dfield<elec_t::coulomb>(                                       //
         r2, xr, yr, zr, dscale, pscale,                                  //
         ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz, pdi, pti, //
         rpole[k][mpl_pme_0], rpole[k][mpl_pme_x], rpole[k][mpl_pme_y],
         rpole[k][mpl_pme_z], rpole[k][mpl_pme_xx], rpole[k][mpl_pme_xy],
         rpole[k][mpl_pme_xz], rpole[k][mpl_pme_yy], rpole[k][mpl_pme_yz],
         rpole[k][mpl_pme_zz], pdamp[k], thole[k], //
         0, pairf);

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
   }
}

// see also subroutine ufield0b in induce.f
void ufield_coulomb(const real (*uind)[3], const real (*uinp)[3],
                    real (*field)[3], real (*fieldp)[3])
{
   device_array::zero(n, field, fieldp);

   const real off = switch_off(switch_mpole);
   const real off2 = off * off;
   const int maxnlst = mlist_unit->maxnlst;
   const auto* mlst = mlist_unit.deviceptr();

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
            PairField pairf;
            pair_ufield<elec_t::coulomb>(                                //
               r2, xr, yr, zr, 1,                                        //
               uindi0, uindi1, uindi2, uinpi0, uinpi1, uinpi2, pdi, pti, //
               uind[k][0], uind[k][1], uind[k][2], uinp[k][0], uinp[k][1],
               uinp[k][2], pdamp[k], thole[k], //
               0, pairf);

            gxi += pairf.fid[0];
            gyi += pairf.fid[1];
            gzi += pairf.fid[2];
            txi += pairf.fip[0];
            tyi += pairf.fip[1];
            tzi += pairf.fip[2];

            atomic_add(pairf.fkd[0], &field[k][0]);
            atomic_add(pairf.fkd[1], &field[k][1]);
            atomic_add(pairf.fkd[2], &field[k][2]);
            atomic_add(pairf.fkp[0], &fieldp[k][0]);
            atomic_add(pairf.fkp[1], &fieldp[k][1]);
            atomic_add(pairf.fkp[2], &fieldp[k][2]);
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
         PairField pairf;
         pair_ufield<elec_t::coulomb>(                                //
            r2, xr, yr, zr, uscale,                                   //
            uindi0, uindi1, uindi2, uinpi0, uinpi1, uinpi2, pdi, pti, //
            uind[k][0], uind[k][1], uind[k][2], uinp[k][0], uinp[k][1],
            uinp[k][2], pdamp[k], thole[k], //
            0, pairf);

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
      }
   }
}
TINKER_NAMESPACE_END
