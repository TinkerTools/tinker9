#include "acc_common.h"
#include "drt_pair_polar.h"
#include "e_polar.h"
#include "gpu_card.h"
#include "md.h"
#include "nblist.h"
#include "pme.h"

TINKER_NAMESPACE_BEGIN
template <int USE>
void epolar_real_tmpl (const real (*uind)[3], const real (*uinp)[3])
{
   constexpr int do_e = USE & calc::energy;
   constexpr int do_a = USE & calc::analyz;
   constexpr int do_g = USE & calc::grad;
   constexpr int do_v = USE & calc::virial;
   sanity_check<USE> ();

   if_constexpr (do_g) device_array::zero (n, ufld, dufld);

   const real off = switch_off (switch_ewald);
   const real off2 = off * off;
   const int maxnlst = mlist_unit->maxnlst;
   const auto* mlst = mlist_unit.deviceptr ();

   auto* nep = ep_handle.ne ()->buffer ();
   auto* ep = ep_handle.e ()->buffer ();
   auto* vir_ep = ep_handle.vir ()->buffer ();
   auto bufsize = ep_handle.buffer_size ();

   const real f = 0.5 * electric / dielec;

   const PMEUnit pu = ppme_unit;
   const real aewald = pu->aewald;

#define POLAR_DPTRS_                                                           \
   x, y, z, gx, gy, gz, box, rpole, thole, pdamp, uind, uinp, nep, ep, vir_ep, \
      ufld, dufld

   MAYBE_UNUSED int GRID_DIM = get_grid_size (BLOCK_DIM);
   #pragma acc parallel num_gangs(GRID_DIM) vector_length(BLOCK_DIM)\
               deviceptr(POLAR_DPTRS_,mlst)
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
      real uix = uind[i][0];
      real uiy = uind[i][1];
      real uiz = uind[i][2];
      real pdi = pdamp[i];
      real pti = thole[i];
      MAYBE_UNUSED real uixp, uiyp, uizp;
      if_constexpr (do_g)
      {
         uixp = uinp[i][0];
         uiyp = uinp[i][1];
         uizp = uinp[i][2];
      }
      real gxi = 0, gyi = 0, gzi = 0;
      real txi = 0, tyi = 0, tzi = 0;
      real du0 = 0, du1 = 0, du2 = 0, du3 = 0, du4 = 0, du5 = 0;

      int nmlsti = mlst->nlst[i];
      int base = i * maxnlst;
      #pragma acc loop vector independent\
                reduction(+:gxi,gyi,gzi,txi,tyi,tzi,du0,du1,du2,du3,du4,du5)
      for (int kk = 0; kk < nmlsti; ++kk) {
         int offset = kk & (bufsize - 1);
         int k = mlst->lst[base + kk];
         real xr = x[k] - xi;
         real yr = y[k] - yi;
         real zr = z[k] - zi;

         image (xr, yr, zr, box);
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
            if_constexpr (do_g)
            {
               ukxp = uinp[k][0];
               ukyp = uinp[k][1];
               ukzp = uinp[k][2];
            }

            MAYBE_UNUSED real e;
            MAYBE_UNUSED PairPolarGrad pgrad;
            pair_polar<USE, elec_t::ewald> ( //
               r2, xr, yr, zr, 1, 1, 1,      //
               ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz, uix, uiy,
               uiz, uixp, uiyp, uizp, pdi, pti, //
               ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukx, uky,
               ukz, ukxp, ukyp, ukzp, pdamp[k], thole[k], //
               f, aewald, e, pgrad);

            if_constexpr (do_a && do_e)
            {
               atomic_add_value (1, nep, offset);
               atomic_add_value (e, ep, offset);
            }

            if_constexpr (do_g)
            {
               gxi += pgrad.frcx;
               gyi += pgrad.frcy;
               gzi += pgrad.frcz;
               atomic_add_value (-pgrad.frcx, gx, k);
               atomic_add_value (-pgrad.frcy, gy, k);
               atomic_add_value (-pgrad.frcz, gz, k);

               txi += pgrad.ufldi[0];
               tyi += pgrad.ufldi[1];
               tzi += pgrad.ufldi[2];
               atomic_add_value (pgrad.ufldk[0], &ufld[k][0]);
               atomic_add_value (pgrad.ufldk[1], &ufld[k][1]);
               atomic_add_value (pgrad.ufldk[2], &ufld[k][2]);

               du0 += pgrad.dufldi[0];
               du1 += pgrad.dufldi[1];
               du2 += pgrad.dufldi[2];
               du3 += pgrad.dufldi[3];
               du4 += pgrad.dufldi[4];
               du5 += pgrad.dufldi[5];
               atomic_add_value (pgrad.dufldk[0], &dufld[k][0]);
               atomic_add_value (pgrad.dufldk[1], &dufld[k][1]);
               atomic_add_value (pgrad.dufldk[2], &dufld[k][2]);
               atomic_add_value (pgrad.dufldk[3], &dufld[k][3]);
               atomic_add_value (pgrad.dufldk[4], &dufld[k][4]);
               atomic_add_value (pgrad.dufldk[5], &dufld[k][5]);

               if_constexpr (do_v)
               {
                  real vxx = -xr * pgrad.frcx;
                  real vxy = -0.5f * (yr * pgrad.frcx + xr * pgrad.frcy);
                  real vxz = -0.5f * (zr * pgrad.frcx + xr * pgrad.frcz);
                  real vyy = -yr * pgrad.frcy;
                  real vyz = -0.5f * (zr * pgrad.frcy + yr * pgrad.frcz);
                  real vzz = -zr * pgrad.frcz;

                  atomic_add_value (vxx, vxy, vxz, vyy, vyz, vzz, vir_ep,
                                    offset);
               }
            }
         }
         // end if use_thole

      } // end for (int kk)

      if_constexpr (do_g)
      {
         atomic_add_value (gxi, gx, i);
         atomic_add_value (gyi, gy, i);
         atomic_add_value (gzi, gz, i);
         atomic_add_value (txi, &ufld[i][0]);
         atomic_add_value (tyi, &ufld[i][1]);
         atomic_add_value (tzi, &ufld[i][2]);
         atomic_add_value (du0, &dufld[i][0]);
         atomic_add_value (du1, &dufld[i][1]);
         atomic_add_value (du2, &dufld[i][2]);
         atomic_add_value (du3, &dufld[i][3]);
         atomic_add_value (du4, &dufld[i][4]);
         atomic_add_value (du5, &dufld[i][5]);
      }
   } // end for (int i)

   #pragma acc parallel deviceptr(POLAR_DPTRS_,dpuexclude_,dpuexclude_scale_)
   #pragma acc loop independent
   for (int ii = 0; ii < ndpuexclude_; ++ii) {
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
      if_constexpr (do_g)
      {
         uixp = uinp[i][0];
         uiyp = uinp[i][1];
         uizp = uinp[i][2];
      }

      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;

      image (xr, yr, zr, box);
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
         if_constexpr (do_g)
         {
            ukxp = uinp[k][0];
            ukyp = uinp[k][1];
            ukzp = uinp[k][2];
         }

         MAYBE_UNUSED real e;
         MAYBE_UNUSED PairPolarGrad pgrad;
         pair_polar<USE, elec_t::coulomb> (         //
            r2, xr, yr, zr, dscale, pscale, uscale, //
            ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz, uix, uiy,
            uiz, uixp, uiyp, uizp, pdi, pti, //
            ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukx, uky,
            ukz, ukxp, ukyp, ukzp, pdamp[k], thole[k], //
            f, 0, e, pgrad);

         if_constexpr (do_a && do_e)
         {
            if (pscale == -1)
               atomic_add_value (-1, nep, offset);
            atomic_add_value (e, ep, offset);
         }

         if_constexpr (do_g)
         {
            atomic_add_value (pgrad.frcx, gx, i);
            atomic_add_value (pgrad.frcy, gy, i);
            atomic_add_value (pgrad.frcz, gz, i);
            atomic_add_value (-pgrad.frcx, gx, k);
            atomic_add_value (-pgrad.frcy, gy, k);
            atomic_add_value (-pgrad.frcz, gz, k);

            atomic_add_value (pgrad.ufldi[0], &ufld[i][0]);
            atomic_add_value (pgrad.ufldi[1], &ufld[i][1]);
            atomic_add_value (pgrad.ufldi[2], &ufld[i][2]);
            atomic_add_value (pgrad.ufldk[0], &ufld[k][0]);
            atomic_add_value (pgrad.ufldk[1], &ufld[k][1]);
            atomic_add_value (pgrad.ufldk[2], &ufld[k][2]);

            atomic_add_value (pgrad.dufldi[0], &dufld[i][0]);
            atomic_add_value (pgrad.dufldi[1], &dufld[i][1]);
            atomic_add_value (pgrad.dufldi[2], &dufld[i][2]);
            atomic_add_value (pgrad.dufldi[3], &dufld[i][3]);
            atomic_add_value (pgrad.dufldi[4], &dufld[i][4]);
            atomic_add_value (pgrad.dufldi[5], &dufld[i][5]);
            atomic_add_value (pgrad.dufldk[0], &dufld[k][0]);
            atomic_add_value (pgrad.dufldk[1], &dufld[k][1]);
            atomic_add_value (pgrad.dufldk[2], &dufld[k][2]);
            atomic_add_value (pgrad.dufldk[3], &dufld[k][3]);
            atomic_add_value (pgrad.dufldk[4], &dufld[k][4]);
            atomic_add_value (pgrad.dufldk[5], &dufld[k][5]);

            if_constexpr (do_v)
            {
               real vxx = -xr * pgrad.frcx;
               real vxy = -0.5f * (yr * pgrad.frcx + xr * pgrad.frcy);
               real vxz = -0.5f * (zr * pgrad.frcx + xr * pgrad.frcz);
               real vyy = -yr * pgrad.frcy;
               real vyz = -0.5f * (zr * pgrad.frcy + yr * pgrad.frcz);
               real vzz = -zr * pgrad.frcz;

               atomic_add_value (vxx, vxy, vxz, vyy, vyz, vzz, vir_ep, offset);
            }
         }
      }
   }

   // torque

   if_constexpr (do_g)
   {
      #pragma acc parallel loop independent\
                deviceptr(rpole,trqx,trqy,trqz,ufld,dufld)
      for (int i = 0; i < n; ++i) {
         real dix = rpole[i][mpl_pme_x];
         real diy = rpole[i][mpl_pme_y];
         real diz = rpole[i][mpl_pme_z];
         real qixx = rpole[i][mpl_pme_xx];
         real qixy = rpole[i][mpl_pme_xy];
         real qixz = rpole[i][mpl_pme_xz];
         real qiyy = rpole[i][mpl_pme_yy];
         real qiyz = rpole[i][mpl_pme_yz];
         real qizz = rpole[i][mpl_pme_zz];

         real tep1 = diz * ufld[i][1] - diy * ufld[i][2] + qixz * dufld[i][1] -
            qixy * dufld[i][3] + 2 * qiyz * (dufld[i][2] - dufld[i][5]) +
            (qizz - qiyy) * dufld[i][4];
         real tep2 = dix * ufld[i][2] - diz * ufld[i][0] - qiyz * dufld[i][1] +
            qixy * dufld[i][4] + 2 * qixz * (dufld[i][5] - dufld[i][0]) +
            (qixx - qizz) * dufld[i][3];
         real tep3 = diy * ufld[i][0] - dix * ufld[i][1] + qiyz * dufld[i][3] -
            qixz * dufld[i][4] + 2 * qixy * (dufld[i][0] - dufld[i][2]) +
            (qiyy - qixx) * dufld[i][1];

         trqx[i] += tep1;
         trqy[i] += tep2;
         trqz[i] += tep3;
      }
   }
}

template <int USE>
void epolar_recip_self_tmpl (const real (*gpu_uind)[3],
                             const real (*gpu_uinp)[3])
{
   constexpr int do_e = USE & calc::energy;
   constexpr int do_a = USE & calc::analyz;
   constexpr int do_g = USE & calc::grad;
   constexpr int do_v = USE & calc::virial;
   sanity_check<USE> ();

   const PMEUnit pu = ppme_unit;
   const auto& st = *pu;
   const int nfft1 = st.nfft1;
   const int nfft2 = st.nfft2;
   const int nfft3 = st.nfft3;
   const real aewald = st.aewald;

   const real f = electric / dielec;

   auto* nep = ep_handle.ne ()->buffer ();
   auto* ep = ep_handle.e ()->buffer ();
   auto bufsize = ep_handle.buffer_size ();

   auto* fphid = fdip_phi1;
   auto* fphip = fdip_phi2;

   cuind_to_fuind (pu, gpu_uind, gpu_uinp, fuind, fuinp);
   if_constexpr (do_e && do_a)
   {
      // if (pairwise .eq. .true.)
      #pragma acc parallel loop independent deviceptr(fuind,fphi,ep)
      for (int i = 0; i < n; ++i) {
         int offset = i & (bufsize - 1);
         real e = 0.5f * f *
            (fuind[i][0] * fphi[i][1] + fuind[i][1] * fphi[i][2] +
             fuind[i][2] * fphi[i][3]);
         atomic_add_value (e, ep, offset);
      }
      // end if
   }
   grid_uind (pu, fuind, fuinp);
   fftfront (pu);
   // TODO: store vs. recompute qfac
   pme_conv0 (pu);
   fftback (pu);
   fphi_uind (pu, fphid, fphip, fphidp);

   // increment the dipole polarization gradient contributions

   // data deriv1  / 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /
   // data deriv2  / 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /
   // data deriv3  / 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /
   constexpr int deriv1[10] = {1, 4, 7, 8, 10, 15, 17, 13, 14, 19};
   constexpr int deriv2[10] = {2, 7, 5, 9, 13, 11, 18, 15, 19, 16};
   constexpr int deriv3[10] = {3, 8, 9, 6, 14, 16, 12, 19, 17, 18};

   #pragma acc parallel loop independent deviceptr(box,gx,gy,gz,\
              fmp,fphi,fuind,fuinp,fphid,fphip,fphidp)
   for (int i = 0; i < n; ++i) {
      real f1 = 0;
      real f2 = 0;
      real f3 = 0;
      #pragma acc loop independent reduction(+:f1,f2,f3)
      for (int k = 0; k < 3; ++k) {
         int j1 = deriv1[k + 1];
         int j2 = deriv2[k + 1];
         int j3 = deriv3[k + 1];
         f1 += (fuind[i][k] + fuinp[i][k]) * fphi[i][j1];
         f2 += (fuind[i][k] + fuinp[i][k]) * fphi[i][j2];
         f3 += (fuind[i][k] + fuinp[i][k]) * fphi[i][j3];
         // if poltyp .eq. 'MUTUAL'
         f1 += fuind[i][k] * fphip[i][j1] + fuinp[i][k] * fphid[i][j1];
         f2 += fuind[i][k] * fphip[i][j2] + fuinp[i][k] * fphid[i][j2];
         f3 += fuind[i][k] * fphip[i][j3] + fuinp[i][k] * fphid[i][j3];
         // end if
      }
      #pragma acc loop independent reduction(+:f1,f2,f3)
      for (int k = 0; k < 10; ++k) {
         f1 += fmp[i][k] * fphidp[i][deriv1[k]];
         f2 += fmp[i][k] * fphidp[i][deriv2[k]];
         f3 += fmp[i][k] * fphidp[i][deriv3[k]];
      }
      if_constexpr (do_g)
      {
         f1 *= 0.5f * nfft1;
         f2 *= 0.5f * nfft2;
         f3 *= 0.5f * nfft3;
         real h1 = box->recip[0][0] * f1 + box->recip[1][0] * f2 +
            box->recip[2][0] * f3;
         real h2 = box->recip[0][1] * f1 + box->recip[1][1] * f2 +
            box->recip[2][1] * f3;
         real h3 = box->recip[0][2] * f1 + box->recip[1][2] * f2 +
            box->recip[2][2] * f3;
         gx[i] += h1 * f;
         gy[i] += h2 * f;
         gz[i] += h3 * f;
      }
   } // end for (int i)

   // set the potential to be the induced dipole average

   // see also subroutine eprecip1 in epolar1.f
   // do i = 1, npole
   //    do j = 1, 10
   //       fphidp(j,i) = 0.5d0 * fphidp(j,i)
   //    end do
   // end do
   // Notice that only 10 * n elements were scaled in the original code.
   device_array::scale (n, 0.5f * f, fphidp);
   fphi_to_cphi (pu, fphidp, cphidp);

   // recip and self torques

   real term = f * REAL_CUBE (aewald) * 4 / 3 / sqrtpi;
   real fterm_term = -2 * f * REAL_CUBE (aewald) / 3 / sqrtpi;
   #pragma acc parallel loop independent\
              deviceptr(ep,nep,trqx,trqy,trqz,\
              rpole,cmp,gpu_uind,gpu_uinp,cphidp)
   for (int i = 0; i < n; ++i) {
      int offset = i & (bufsize - 1);
      real dix = rpole[i][mpl_pme_x];
      real diy = rpole[i][mpl_pme_y];
      real diz = rpole[i][mpl_pme_z];
      real uix = 0.5f * (gpu_uind[i][0] + gpu_uinp[i][0]);
      real uiy = 0.5f * (gpu_uind[i][1] + gpu_uinp[i][1]);
      real uiz = 0.5f * (gpu_uind[i][2] + gpu_uinp[i][2]);

      if_constexpr (do_g)
      {
         real tep1 = cmp[i][3] * cphidp[i][2] - cmp[i][2] * cphidp[i][3] +
            2 * (cmp[i][6] - cmp[i][5]) * cphidp[i][9] +
            cmp[i][8] * cphidp[i][7] + cmp[i][9] * cphidp[i][5] -
            cmp[i][7] * cphidp[i][8] - cmp[i][9] * cphidp[i][6];
         real tep2 = cmp[i][1] * cphidp[i][3] - cmp[i][3] * cphidp[i][1] +
            2 * (cmp[i][4] - cmp[i][6]) * cphidp[i][8] +
            cmp[i][7] * cphidp[i][9] + cmp[i][8] * cphidp[i][6] -
            cmp[i][8] * cphidp[i][4] - cmp[i][9] * cphidp[i][7];
         real tep3 = cmp[i][2] * cphidp[i][1] - cmp[i][1] * cphidp[i][2] +
            2 * (cmp[i][5] - cmp[i][4]) * cphidp[i][7] +
            cmp[i][7] * cphidp[i][4] + cmp[i][9] * cphidp[i][8] -
            cmp[i][7] * cphidp[i][5] - cmp[i][8] * cphidp[i][9];

         // self term

         tep1 += term * (diy * uiz - diz * uiy);
         tep2 += term * (diz * uix - dix * uiz);
         tep3 += term * (dix * uiy - diy * uix);

         trqx[i] += tep1;
         trqy[i] += tep2;
         trqz[i] += tep3;
      }

      if_constexpr (do_e && do_a)
      {
         // if (pairwise .eq. .true.)
         uix = gpu_uind[i][0];
         uiy = gpu_uind[i][1];
         uiz = gpu_uind[i][2];
         real uii = dix * uix + diy * uiy + diz * uiz;
         atomic_add_value (fterm_term * uii, ep, offset);
         atomic_add_value (1, nep, offset);
         // end if
      }
   }

   // recip virial

   if_constexpr (do_v)
   {

      auto* vir_ep = ep_handle.vir ()->buffer ();
      auto* vir_m = vir_m_handle->buffer ();
      auto vir_m_len = vir_m_handle->size ();
      assert (bufsize >= vir_m_len);

      #pragma acc parallel loop independent deviceptr(vir_ep,vir_m)
      for (int i = 0; i < vir_m_len * VirialBuffer::NS; ++i) {
         vir_ep[i] -= vir_m[i];
      }

      device_array::scale (n, f, cphi, fphid, fphip);

      #pragma acc parallel loop independent\
                deviceptr(vir_ep,box,cmp,\
                gpu_uind,gpu_uinp,fphid,fphip,cphi,cphidp)
      for (int i = 0; i < n; ++i) {
         real cphid[4], cphip[4];
         real ftc[3][3];

         // frac_to_cart

         ftc[0][0] = nfft1 * box->recip[0][0];
         ftc[1][0] = nfft2 * box->recip[1][0];
         ftc[2][0] = nfft3 * box->recip[2][0];
         ftc[0][1] = nfft1 * box->recip[0][1];
         ftc[1][1] = nfft2 * box->recip[1][1];
         ftc[2][1] = nfft3 * box->recip[2][1];
         ftc[0][2] = nfft1 * box->recip[0][2];
         ftc[1][2] = nfft2 * box->recip[1][2];
         ftc[2][2] = nfft3 * box->recip[2][2];

         #pragma acc loop independent
         for (int j = 0; j < 3; ++j) {
            cphid[j + 1] = 0;
            cphip[j + 1] = 0;
            #pragma acc loop seq
            for (int k = 0; k < 3; ++k) {
               cphid[j + 1] += ftc[k][j] * fphid[i][k + 1];
               cphip[j + 1] += ftc[k][j] * fphip[i][k + 1];
            }
         }

         real vxx = 0;
         real vyy = 0;
         real vzz = 0;
         real vxy = 0;
         real vxz = 0;
         real vyz = 0;

         vxx = vxx - cmp[i][1] * cphidp[i][1] -
            0.5f * ((gpu_uind[i][0] + gpu_uinp[i][0]) * cphi[i][1]);
         vxy = vxy -
            0.5f * (cphidp[i][1] * cmp[i][2] + cphidp[i][2] * cmp[i][1]) -
            0.25f *
               ((gpu_uind[i][1] + gpu_uinp[i][1]) * cphi[i][1] +
                (gpu_uind[i][0] + gpu_uinp[i][0]) * cphi[i][2]);
         vxz = vxz -
            0.5f * (cphidp[i][1] * cmp[i][3] + cphidp[i][3] * cmp[i][1]) -
            0.25f *
               ((gpu_uind[i][2] + gpu_uinp[i][2]) * cphi[i][1] +
                (gpu_uind[i][0] + gpu_uinp[i][0]) * cphi[i][3]);
         vyy = vyy - cmp[i][2] * cphidp[i][2] -
            0.5f * ((gpu_uind[i][1] + gpu_uinp[i][1]) * cphi[i][2]);
         vyz = vyz -
            0.5f * (cphidp[i][2] * cmp[i][3] + cphidp[i][3] * cmp[i][2]) -
            0.25f *
               ((gpu_uind[i][2] + gpu_uinp[i][2]) * cphi[i][2] +
                (gpu_uind[i][1] + gpu_uinp[i][1]) * cphi[i][3]);
         vzz = vzz - cmp[i][3] * cphidp[i][3] -
            0.5f * ((gpu_uind[i][2] + gpu_uinp[i][2]) * cphi[i][3]);
         vxx = vxx - 2 * cmp[i][4] * cphidp[i][4] - cmp[i][7] * cphidp[i][7] -
            cmp[i][8] * cphidp[i][8];
         vxy = vxy - (cmp[i][4] + cmp[i][5]) * cphidp[i][7] -
            0.5f *
               (cmp[i][7] * (cphidp[i][5] + cphidp[i][4]) +
                cmp[i][8] * cphidp[i][9] + cmp[i][9] * cphidp[i][8]);
         vxz = vxz - (cmp[i][4] + cmp[i][6]) * cphidp[i][8] -
            0.5f *
               (cmp[i][8] * (cphidp[i][4] + cphidp[i][6]) +
                cmp[i][7] * cphidp[i][9] + cmp[i][9] * cphidp[i][7]);
         vyy = vyy - 2 * cmp[i][5] * cphidp[i][5] - cmp[i][7] * cphidp[i][7] -
            cmp[i][9] * cphidp[i][9];
         vyz = vyz - (cmp[i][5] + cmp[i][6]) * cphidp[i][9] -
            0.5f *
               (cmp[i][9] * (cphidp[i][5] + cphidp[i][6]) +
                cmp[i][7] * cphidp[i][8] + cmp[i][8] * cphidp[i][7]);
         vzz = vzz - 2 * cmp[i][6] * cphidp[i][6] - cmp[i][8] * cphidp[i][8] -
            cmp[i][9] * cphidp[i][9];

         // if (poltyp .eq. 'MUTUAL')
         vxx = vxx -
            0.5f * (cphid[1] * gpu_uinp[i][0] + cphip[1] * gpu_uind[i][0]);
         vxy = vxy -
            0.25f *
               (cphid[1] * gpu_uinp[i][1] + cphip[1] * gpu_uind[i][1] +
                cphid[2] * gpu_uinp[i][0] + cphip[2] * gpu_uind[i][0]);
         vxz = vxz -
            0.25f *
               (cphid[1] * gpu_uinp[i][2] + cphip[1] * gpu_uind[i][2] +
                cphid[3] * gpu_uinp[i][0] + cphip[3] * gpu_uind[i][0]);
         vyy = vyy -
            0.5f * (cphid[2] * gpu_uinp[i][1] + cphip[2] * gpu_uind[i][1]);
         vyz = vyz -
            0.25f *
               (cphid[2] * gpu_uinp[i][2] + cphip[2] * gpu_uind[i][2] +
                cphid[3] * gpu_uinp[i][1] + cphip[3] * gpu_uind[i][1]);
         vzz = vzz -
            0.5f * (cphid[3] * gpu_uinp[i][2] + cphip[3] * gpu_uind[i][2]);
         // end if

         atomic_add_value (vxx, vxy, vxz, vyy, vyz, vzz, vir_ep,
                           i & (bufsize - 1));
      }

      // qgrip: pvu_qgrid
      const PMEUnit pvu = pvpme_unit;
      #pragma acc parallel loop independent deviceptr(cmp,gpu_uinp)
      for (int i = 0; i < n; ++i) {
         cmp[i][1] += gpu_uinp[i][0];
         cmp[i][2] += gpu_uinp[i][1];
         cmp[i][3] += gpu_uinp[i][2];
      }
      cmp_to_fmp (pvu, cmp, fmp);
      grid_mpole (pvu, fmp);
      fftfront (pvu);

      // qgrid: pu_qgrid
      #pragma acc parallel loop independent deviceptr(cmp,gpu_uind,gpu_uinp)
      for (int i = 0; i < n; ++i) {
         cmp[i][1] += (gpu_uind[i][0] - gpu_uinp[i][0]);
         cmp[i][2] += (gpu_uind[i][1] - gpu_uinp[i][1]);
         cmp[i][3] += (gpu_uind[i][2] - gpu_uinp[i][2]);
      }
      cmp_to_fmp (pu, cmp, fmp);
      grid_mpole (pu, fmp);
      fftfront (pu);

      const auto* d = pu.deviceptr ();
      const auto* p = pvu.deviceptr ();
      const int nff = nfft1 * nfft2;
      const int ntot = nfft1 * nfft2 * nfft3;
      real pterm = REAL_SQ (pi / aewald);

      #pragma acc parallel loop independent\
                deviceptr(box,d,p,vir_ep)
      for (int i = 1; i < ntot; ++i) {
         const real volterm = pi * box->volbox;

         int k3 = i / nff;
         int j = i - k3 * nff;
         int k2 = j / nfft1;
         int k1 = j - k2 * nfft1;

         int r1 = (k1 < (nfft1 + 1) / 2) ? k1 : (k1 - nfft1);
         int r2 = (k2 < (nfft2 + 1) / 2) ? k2 : (k2 - nfft2);
         int r3 = (k3 < (nfft3 + 1) / 2) ? k3 : (k3 - nfft3);

         real h1 = box->recip[0][0] * r1 + box->recip[1][0] * r2 +
            box->recip[2][0] * r3;
         real h2 = box->recip[0][1] * r1 + box->recip[1][1] * r2 +
            box->recip[2][1] * r3;
         real h3 = box->recip[0][2] * r1 + box->recip[1][2] * r2 +
            box->recip[2][2] * r3;
         real hsq = h1 * h1 + h2 * h2 + h3 * h3;
         real term = -pterm * hsq;
         real expterm = 0;
         if (term > -50) {
            // TODO: if .not. use_bounds; if octahedron; 2/hsq
            real denom =
               volterm * hsq * d->bsmod1[k1] * d->bsmod1[k2] * d->bsmod1[k3];
            expterm = REAL_EXP (term) / denom;

            real struc2 = d->qgrid[2 * i] * p->qgrid[2 * i] +
               d->qgrid[2 * i + 1] * p->qgrid[2 * i + 1];
            real eterm = 0.5f * f * expterm * struc2;
            real vterm = (2 / hsq) * (1 - term) * eterm;

            real vxx = (h1 * h1 * vterm - eterm);
            real vxy = h1 * h2 * vterm;
            real vxz = h1 * h3 * vterm;
            real vyy = (h2 * h2 * vterm - eterm);
            real vyz = h2 * h3 * vterm;
            real vzz = (h3 * h3 * vterm - eterm);

            atomic_add_value (vxx, vxy, vxz, vyy, vyz, vzz, vir_ep,
                              i & (bufsize - 1));
         }
      }
   }
}

template <int USE>
void epolar_ewald_tmpl (const real (*uind)[3], const real (*uinp)[3])
{
   constexpr int do_e = USE & calc::energy;
   constexpr int do_a = USE & calc::analyz;
   constexpr int do_g = USE & calc::grad;
   sanity_check<USE> ();

   if_constexpr (do_e && !do_a) epolar0_dotprod (uind, udirp);
   static_assert (do_g || do_a,
                  "Do not use this template for the energy-only version.");

   epolar_real_tmpl<USE> (uind, uinp);

   epolar_recip_self_tmpl<USE> (uind, uinp);
}

void epolar_ewald (int vers)
{
   if (vers == calc::v0) {
      induce (uind, uinp);
      epolar0_dotprod (uind, udirp);
   } else if (vers == calc::v1) {
      induce (uind, uinp);
      epolar_ewald_tmpl<calc::v1> (uind, uinp);
   } else if (vers == calc::v3) {
      induce (uind, uinp);
      epolar_ewald_tmpl<calc::v3> (uind, uinp);
   } else if (vers == calc::v4) {
      induce (uind, uinp);
      epolar_ewald_tmpl<calc::v4> (uind, uinp);
   } else if (vers == calc::v5) {
      induce (uind, uinp);
      epolar_ewald_tmpl<calc::v5> (uind, uinp);
   } else if (vers == calc::v6) {
      induce (uind, uinp);
      epolar_ewald_tmpl<calc::v6> (uind, uinp);
   }
}
TINKER_NAMESPACE_END
