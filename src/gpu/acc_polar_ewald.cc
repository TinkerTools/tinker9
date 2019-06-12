#include "acc_e.h"
#include "gpu/decl_pme.h"
#include "gpu/e_polar.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
template <int USE>
void epolar_real_tmpl(const real (*gpu_uind)[3], const real (*gpu_uinp)[3]) {}

template <int USE>
void epolar_recip_self_tmpl(const real (*gpu_uind)[3],
                            const real (*gpu_uinp)[3]) {
  const int pu = ppme_unit;
  const auto& st = pme_obj(pu);
  const int nfft1 = st.nfft1;
  const int nfft2 = st.nfft2;
  const int nfft3 = st.nfft3;
  const real aewald = st.aewald;

  real(*fphid)[10] = fdip_phi1;
  real(*fphip)[10] = fdip_phi2;

  cuind_to_fuind(pu, gpu_uind, gpu_uinp, fuind, fuinp);
  grid_uind(pu, fuind, fuinp);
  fftfront(pu);
  // TODO: store vs. recompute qfac
  pme_conv0(pu);
  fftback(pu);
  fphi_uind(pu, fphid, fphip, fphidp);

  const real f = chgpot::electric / chgpot::dielec;

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
    f1 *= 0.5f * nfft1;
    f2 *= 0.5f * nfft2;
    f3 *= 0.5f * nfft3;
    real h1 =
        box->recip[0][0] * f1 + box->recip[1][0] * f2 + box->recip[2][0] * f3;
    real h2 =
        box->recip[0][1] * f1 + box->recip[1][1] * f2 + box->recip[2][1] * f3;
    real h3 =
        box->recip[0][2] * f1 + box->recip[1][2] * f2 + box->recip[2][2] * f3;
    gx[i] += h1 * f;
    gy[i] += h2 * f;
    gz[i] += h3 * f;
  } // end for (int i)

  // set the potential to be the induced dipole average

  // see also subroutine eprecip1 in epolar1.f
  // do i = 1, npole
  //    do j = 1, 10
  //       fphidp(j,i) = 0.5d0 * fphidp(j,i)
  //    end do
  // end do
  // Notice that only 10 * n elements were scaled in the original code.
  scale_data(&fphidp[0][0], 0.5f * f, 20 * n);
  fphi_to_cphi(pu, fphidp, cphidp);

  // recip and self torques

  real term = f * REAL_CUBE(aewald) * 4 / 3 / sqrtpi;
  real cmp[10];
  #pragma acc parallel loop independent deviceptr(trqx,trqy,trqz,\
              rpole,gpu_uind,gpu_uinp,cphidp) private(cmp[0:10])
  for (int i = 0; i < n; ++i) {
    cmp[0] = rpole[i][mpl_pme_0];
    cmp[1] = rpole[i][mpl_pme_x];
    cmp[2] = rpole[i][mpl_pme_y];
    cmp[3] = rpole[i][mpl_pme_z];
    cmp[4] = rpole[i][mpl_pme_xx];
    cmp[5] = rpole[i][mpl_pme_yy];
    cmp[6] = rpole[i][mpl_pme_zz];
    cmp[7] = 2 * rpole[i][mpl_pme_xy];
    cmp[8] = 2 * rpole[i][mpl_pme_xz];
    cmp[9] = 2 * rpole[i][mpl_pme_yz];
    real tep1 = cmp[3] * cphidp[i][2] - cmp[2] * cphidp[i][3] +
        2 * (cmp[6] - cmp[5]) * cphidp[i][9] + cmp[8] * cphidp[i][7] +
        cmp[9] * cphidp[i][5] - cmp[7] * cphidp[i][8] - cmp[9] * cphidp[i][6];
    real tep2 = cmp[1] * cphidp[i][3] - cmp[3] * cphidp[i][1] +
        2 * (cmp[4] - cmp[6]) * cphidp[i][8] + cmp[7] * cphidp[i][9] +
        cmp[8] * cphidp[i][6] - cmp[8] * cphidp[i][4] - cmp[9] * cphidp[i][7];
    real tep3 = cmp[2] * cphidp[i][1] - cmp[1] * cphidp[i][2] +
        2 * (cmp[5] - cmp[4]) * cphidp[i][7] + cmp[7] * cphidp[i][4] +
        cmp[9] * cphidp[i][8] - cmp[7] * cphidp[i][5] - cmp[8] * cphidp[i][9];

    real dix = rpole[i][mpl_pme_x];
    real diy = rpole[i][mpl_pme_y];
    real diz = rpole[i][mpl_pme_z];
    real uix = 0.5f * (gpu_uind[i][0] + gpu_uinp[i][0]);
    real uiy = 0.5f * (gpu_uind[i][1] + gpu_uinp[i][1]);
    real uiz = 0.5f * (gpu_uind[i][2] + gpu_uinp[i][2]);
    tep1 += term * (diy * uiz - diz * uiy);
    tep2 += term * (diz * uix - dix * uiz);
    tep3 += term * (dix * uiy - diy * uix);

    trqx[i] += tep1;
    trqy[i] += tep2;
    trqz[i] += tep3;
  }
}

template <int USE>
void epolar_ewald_tmpl(const real (*gpu_uind)[3], const real (*gpu_uinp)[3]) {
  constexpr int do_e = USE & use_energy;
  constexpr int do_a = USE & use_analyz;
  constexpr int do_g = USE & use_grad;
  constexpr int do_v = USE & use_virial;
  static_assert(do_v ? do_g : true, "");
  static_assert(do_a ? do_e : true, "");

  if_constexpr(do_e || do_a || do_v) {
    #pragma acc serial deviceptr(ep,nep,vir_ep)
    {
      if_constexpr(do_e) { *ep = 0; }
      if_constexpr(do_a) { *nep = 0; }
      if_constexpr(do_v) {
        for (int i = 0; i < 9; ++i) {
          vir_ep[i] = 0;
        }
      }
    }
  }

  if_constexpr(do_e && !do_a) epolar0_dotprod(gpu_uind, udirp);
  static_assert(do_g || do_a,
                "Do not use this template for the energy-only version.");

  epolar_real_tmpl<USE>(gpu_uind, gpu_uinp);

  epolar_recip_self_tmpl<USE>(gpu_uind, gpu_uinp);
}
}
TINKER_NAMESPACE_END

extern "C" {
m_tinker_using_namespace;
void tinker_gpu_epolar_ewald0() {
  gpu::induce(&gpu::uind[0][0], &gpu::uinp[0][0]);
  gpu::epolar0_dotprod(gpu::uind, gpu::udirp);
}
void tinker_gpu_epolar_ewald1() {
  gpu::induce(&gpu::uind[0][0], &gpu::uinp[0][0]);
  gpu::epolar_ewald_tmpl<gpu::v1>(gpu::uind, gpu::uinp);
}
void tinker_gpu_epolar_ewald3() {
  gpu::induce(&gpu::uind[0][0], &gpu::uinp[0][0]);
  gpu::epolar_ewald_tmpl<gpu::v3>(gpu::uind, gpu::uinp);
}
void tinker_gpu_epolar_ewald4() {
  gpu::induce(&gpu::uind[0][0], &gpu::uinp[0][0]);
  gpu::epolar_ewald_tmpl<gpu::v4>(gpu::uind, gpu::uinp);
}
void tinker_gpu_epolar_ewald5() {
  gpu::induce(&gpu::uind[0][0], &gpu::uinp[0][0]);
  gpu::epolar_ewald_tmpl<gpu::v5>(gpu::uind, gpu::uinp);
}
void tinker_gpu_epolar_ewald6() {
  gpu::induce(&gpu::uind[0][0], &gpu::uinp[0][0]);
  gpu::epolar_ewald_tmpl<gpu::v6>(gpu::uind, gpu::uinp);
}
}
