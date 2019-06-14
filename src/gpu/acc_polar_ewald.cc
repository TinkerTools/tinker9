#include "acc_e.h"
#include "gpu/decl_pme.h"
#include "gpu/e_polar.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
template <int USE>
void epolar_real_tmpl(const real (*gpu_uind)[3], const real (*gpu_uinp)[3]) {
  constexpr int do_e = USE & use_energy;
  constexpr int do_a = USE & use_analyz;
  constexpr int do_g = USE & use_grad;
  constexpr int do_v = USE & use_virial;
  sanity_check<USE>();

  if_constexpr(do_g) {
    zero_data(&ufld[0][0], 3 * n);
    zero_data(&dufld[0][0], 6 * n);
  }

  const real(*uind)[3] = reinterpret_cast<const real(*)[3]>(gpu_uind);
  const real(*uinp)[3] = reinterpret_cast<const real(*)[3]>(gpu_uinp);

  const real off = ewald_switch_off;
  const real off2 = off * off;
  const int maxnlst = mlist_obj_.maxnlst;

  const real p2scale = polpot::p2scale;
  const real p3scale = polpot::p3scale;
  const real p4scale = polpot::p4scale;
  const real p5scale = polpot::p5scale;

  const real p2iscale = polpot::p2iscale;
  const real p3iscale = polpot::p3iscale;
  const real p4iscale = polpot::p4iscale;
  const real p5iscale = polpot::p5iscale;

  const real d1scale = polpot::d1scale;
  const real d2scale = polpot::d2scale;
  const real d3scale = polpot::d3scale;
  const real d4scale = polpot::d4scale;

  const real u1scale = polpot::u1scale;
  const real u2scale = polpot::u2scale;
  const real u3scale = polpot::u3scale;
  const real u4scale = polpot::u4scale;

  static std::vector<real> pscalebuf;
  static std::vector<real> dscalebuf;
  static std::vector<real> uscalebuf;
  pscalebuf.resize(n, 1);
  dscalebuf.resize(n, 1);
  uscalebuf.resize(n, 1);
  real* pscale = pscalebuf.data();
  real* dscale = dscalebuf.data();
  real* uscale = uscalebuf.data();

  const real f = 0.5 * chgpot::electric / chgpot::dielec;

  const real aewald = pme_obj(ppme_unit).aewald;
  real bn[5];

  #pragma acc parallel loop independent\
              deviceptr(x,y,z,box,couple,polargroup,mlst,\
              rpole,thole,pdamp,uind,uinp,\
              ep,nep,vir_ep,ufld,dufld)\
              firstprivate(pscale[0:n],dscale[0:n],uscale[0:n])
  for (int i = 0; i < n; ++i) {

    // set exclusion coefficients for connected atoms

    const int n12i = couple->n12[i];
    const int n13i = couple->n13[i];
    const int n14i = couple->n14[i];
    const int n15i = couple->n15[i];

    const int np11i = polargroup->np11[i];
    const int np12i = polargroup->np12[i];
    const int np13i = polargroup->np13[i];
    const int np14i = polargroup->np14[i];

    #pragma acc loop independent
    for (int j = 0; j < n12i; ++j) {
      int ij = couple->i12[i][j];
      pscale[ij] = p2scale;
      #pragma acc loop independent
      for (int k = 0; k < np11i; ++k)
        if (ij == polargroup->ip11[i][k])
          pscale[ij] = p2iscale;
    }
    #pragma acc loop independent
    for (int j = 0; j < n13i; ++j) {
      int ij = couple->i13[i][j];
      pscale[ij] = p3scale;
      #pragma acc loop independent
      for (int k = 0; k < np11i; ++k)
        if (ij == polargroup->ip11[i][k])
          pscale[ij] = p3iscale;
    }
    #pragma acc loop independent
    for (int j = 0; j < n14i; ++j) {
      int ij = couple->i14[i][j];
      pscale[ij] = p4scale;
      #pragma acc loop independent
      for (int k = 0; k < np11i; ++k)
        if (ij == polargroup->ip11[i][k])
          pscale[ij] = p4iscale;
    }
    #pragma acc loop independent
    for (int j = 0; j < n15i; ++j) {
      int ij = couple->i15[i][j];
      pscale[ij] = p5scale;
      #pragma acc loop independent
      for (int k = 0; k < np11i; ++k)
        if (ij == polargroup->ip11[i][k])
          pscale[ij] = p5iscale;
    }

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
    if_constexpr(do_g) {
      uixp = uinp[i][0];
      uiyp = uinp[i][1];
      uizp = uinp[i][2];
    }

    int nmlsti = mlst->nlst[i];
    int base = i * maxnlst;
    #pragma acc loop independent private(bn[0:5])
    for (int kk = 0; kk < nmlsti; ++kk) {
      int k = mlst->lst[base + kk];
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
        if_constexpr(do_g) {
          ukxp = uinp[k][0];
          ukyp = uinp[k][1];
          ukzp = uinp[k][2];
        }

        real dir = dix * xr + diy * yr + diz * zr;
        real qix = qixx * xr + qixy * yr + qixz * zr;
        real qiy = qixy * xr + qiyy * yr + qiyz * zr;
        real qiz = qixz * xr + qiyz * yr + qizz * zr;
        real qir = qix * xr + qiy * yr + qiz * zr;
        real dkr = dkx * xr + dky * yr + dkz * zr;
        real qkx = qkxx * xr + qkxy * yr + qkxz * zr;
        real qky = qkxy * xr + qkyy * yr + qkyz * zr;
        real qkz = qkxz * xr + qkyz * yr + qkzz * zr;
        real qkr = qkx * xr + qky * yr + qkz * zr;
        real diu = dix * ukx + diy * uky + diz * ukz;
        real qiu = qix * ukx + qiy * uky + qiz * ukz;
        real uir = uix * xr + uiy * yr + uiz * zr;
        real dku = dkx * uix + dky * uiy + dkz * uiz;
        real qku = qkx * uix + qky * uiy + qkz * uiz;
        real ukr = ukx * xr + uky * yr + ukz * zr;

        real r = REAL_SQRT(r2);
        real invr1 = REAL_RECIP(r);
        real rr2 = invr1 * invr1;

        real rr1 = f * invr1;
        real rr3 = rr1 * rr2;
        real rr5 = 3 * rr3 * rr2;
        real rr7 = 5 * rr5 * rr2;
        MAYBE_UNUSED real rr9;

        real ralpha = aewald * r;
        bn[0] = REAL_ERFC(ralpha) * invr1;
        real alsq2 = 2 * REAL_SQ(aewald);
        real alsq2n = (aewald > 0 ? REAL_RECIP(sqrtpi * aewald) : 0);
        real exp2a = REAL_EXP(-REAL_SQ(ralpha));
        if_constexpr(!do_g) {
          #pragma acc loop seq
          for (int j = 1; j < 4; ++j) {
            alsq2n *= alsq2;
            bn[j] = ((j + j - 1) * bn[j - 1] + alsq2n * exp2a) * rr2;
          }
        }
        else {
          #pragma acc loop seq
          for (int j = 1; j < 3; ++j) {
            alsq2n *= alsq2;
            bn[j] = ((j + j - 1) * bn[j - 1] + alsq2n * exp2a) * rr2;
          }
        }
        bn[0] *= f;
        bn[1] *= f;
        bn[2] *= f;
        bn[3] *= f;
        if_constexpr(do_g) {
          bn[4] *= f;
          rr9 = 7 * rr7 * rr2;
        }

        real sc3 = 1;
        real sc5 = 1;
        real sc7 = 1;
        MAYBE_UNUSED real rc31, rc32, rc33, rc51, rc52, rc53, rc71, rc72, rc73;
        if_constexpr(do_g) {
          rc31 = 0;
          rc32 = 0;
          rc33 = 0;
          rc51 = 0;
          rc52 = 0;
          rc53 = 0;
          rc71 = 0;
          rc72 = 0;
          rc73 = 0;
        }

        // if use_thole
        real damp = pdi * pdamp[k];
        if (damp != 0) {
          real pgamma = REAL_MIN(pti, thole[k]);
          damp = -pgamma * REAL_CUBE(r / damp);
          if (damp > -50) {
            real expdamp = REAL_EXP(damp);
            sc3 = 1 - expdamp;
            sc5 = 1 - expdamp * (1 - damp);
            sc7 = 1 - expdamp * (1 - damp + (real)0.6 * REAL_SQ(damp));
            if_constexpr(do_g) {
              real temp3 = -3 * damp * expdamp * rr2;
              real temp5 = -damp;
              real temp7 = ((real)-0.2) - ((real)0.6) * damp;
              rc31 = xr * temp3;
              rc32 = yr * temp3;
              rc33 = zr * temp3;
              rc51 = rc31 * temp5;
              rc52 = rc32 * temp5;
              rc53 = rc33 * temp5;
              rc71 = rc51 * temp7;
              rc72 = rc52 * temp7;
              rc73 = rc53 * temp7;
            }
          }
        }

        if_constexpr(do_e && do_a) {
          // if (pairwise .eq. .true.)
          real scalek = pscale[k];
          real sr3 = scalek * sc3 * rr3;
          real sr5 = scalek * sc5 * rr5;
          real sr7 = scalek * sc7 * rr7;
          real term1 = ck * uir - ci * ukr + diu + dku;
          real term2 = 2 * (qiu - qku) - uir * dkr - dir * ukr;
          real term3 = uir * qkr - ukr * qir;
          real efull = term1 * sr3 + term2 * sr5 + term3 * sr7;
          sr3 += (bn[1] - rr3);
          sr5 += (bn[2] - rr5);
          sr7 += (bn[3] - rr7);
          real e = term1 * sr3 + term2 * sr5 + term3 * sr7;
          #pragma acc atomic update
          *ep += e;
          if (efull != 0) {
            #pragma acc atomic update
            *nep += 1;
          }
          // end if
        }
        // end if use_thole
      }
    } // end for (int kk)

    #pragma acc loop independent
    for (int j = 0; j < np11i; ++j) {
      dscale[polargroup->ip11[i][j]] = d1scale;
      uscale[polargroup->ip11[i][j]] = u1scale;
    }
    #pragma acc loop independent
    for (int j = 0; j < np12i; ++j) {
      dscale[polargroup->ip12[i][j]] = d2scale;
      uscale[polargroup->ip12[i][j]] = u2scale;
    }
    #pragma acc loop independent
    for (int j = 0; j < np13i; ++j) {
      dscale[polargroup->ip13[i][j]] = d3scale;
      uscale[polargroup->ip13[i][j]] = u3scale;
    }
    #pragma acc loop independent
    for (int j = 0; j < np14i; ++j) {
      dscale[polargroup->ip14[i][j]] = d4scale;
      uscale[polargroup->ip14[i][j]] = u4scale;
    }

    // reset exclusion coefficients for connected atoms

    #pragma acc loop independent
    for (int j = 0; j < n12i; ++j)
      pscale[couple->i12[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < n13i; ++j)
      pscale[couple->i13[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < n14i; ++j)
      pscale[couple->i14[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < n15i; ++j)
      pscale[couple->i15[i][j]] = 1;

    #pragma acc loop independent
    for (int j = 0; j < np11i; ++j) {
      dscale[polargroup->ip11[i][j]] = 1;
      uscale[polargroup->ip11[i][j]] = 1;
    }
    #pragma acc loop independent
    for (int j = 0; j < np12i; ++j) {
      dscale[polargroup->ip12[i][j]] = 1;
      uscale[polargroup->ip12[i][j]] = 1;
    }
    #pragma acc loop independent
    for (int j = 0; j < np13i; ++j) {
      dscale[polargroup->ip13[i][j]] = 1;
      uscale[polargroup->ip13[i][j]] = 1;
    }
    #pragma acc loop independent
    for (int j = 0; j < np14i; ++j) {
      dscale[polargroup->ip14[i][j]] = 1;
      uscale[polargroup->ip14[i][j]] = 1;
    }
  } // end for (int i)
}

template <int USE>
void epolar_recip_self_tmpl(const real (*gpu_uind)[3],
                            const real (*gpu_uinp)[3]) {
  constexpr int do_e = USE & use_energy;
  constexpr int do_a = USE & use_analyz;
  constexpr int do_g = USE & use_grad;
  constexpr int do_v = USE & use_virial;
  sanity_check<USE>();

  const int pu = ppme_unit;
  const auto& st = pme_obj(pu);
  const int nfft1 = st.nfft1;
  const int nfft2 = st.nfft2;
  const int nfft3 = st.nfft3;
  const real aewald = st.aewald;

  const real f = chgpot::electric / chgpot::dielec;

  real(*fphid)[10] = fdip_phi1;
  real(*fphip)[10] = fdip_phi2;

  cuind_to_fuind(pu, gpu_uind, gpu_uinp, fuind, fuinp);
  if (do_e && do_a) {
    // if (pairwise .eq. .true.)
    #pragma acc parallel loop independent deviceptr(fuind,fphi,ep)
    for (int i = 0; i < n; ++i) {
      real e = 0.5f * f *
          (fuind[i][0] * fphi[i][1] + fuind[i][1] * fphi[i][2] +
           fuind[i][2] * fphi[i][3]);
      #pragma acc atomic update
      *ep += e;
    }
    // end if
  }
  grid_uind(pu, fuind, fuinp);
  fftfront(pu);
  // TODO: store vs. recompute qfac
  pme_conv0(pu);
  fftback(pu);
  fphi_uind(pu, fphid, fphip, fphidp);

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
  real fterm_term = -2 * f * REAL_CUBE(aewald) / 3 / sqrtpi;
  real cmp[10];
  #pragma acc parallel loop independent\
              deviceptr(ep,nep,trqx,trqy,trqz,\
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

    // self term

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

    if_constexpr(do_e && do_a) {
      // if (pairwise .eq. .true.)
      uix = gpu_uind[i][0];
      uiy = gpu_uind[i][1];
      uiz = gpu_uind[i][2];
      real uii = dix * uix + diy * uiy + diz * uiz;
      #pragma acc atomic update
      *ep += fterm_term * uii;
      #pragma acc atomic update
      *nep += 1;
      // end if
    }
  }
}

template <int USE>
void epolar_ewald_tmpl(const real (*gpu_uind)[3], const real (*gpu_uinp)[3]) {
  constexpr int do_e = USE & use_energy;
  constexpr int do_a = USE & use_analyz;
  constexpr int do_g = USE & use_grad;
  constexpr int do_v = USE & use_virial;
  sanity_check<USE>();

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
