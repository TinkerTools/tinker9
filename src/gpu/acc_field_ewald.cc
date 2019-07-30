#include "gpu/acc.h"
#include "gpu/decl_mdstate.h"
#include "gpu/decl_pme.h"
#include "gpu/e_mpole.h"
#include "gpu/e_polar.h"
#include "util_potent.h"

TINKER_NAMESPACE_BEGIN

namespace gpu {
// see also subroutine udirect1 in induce.f
void dfield_ewald_recip_self(real* gpu_field) {
  const int pu = ppme_unit;
  const real aewald = pme_obj(pu).aewald;
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

  real(*field)[3] = reinterpret_cast<real(*)[3]>(gpu_field);

  #pragma acc parallel loop independent deviceptr(field,cphi,rpole)
  for (int i = 0; i < n; ++i) {
    real dix = rpole[i][mpl_pme_x];
    real diy = rpole[i][mpl_pme_y];
    real diz = rpole[i][mpl_pme_z];

    // increment the field at each multipole site
    // get the self-energy portion of the permanent field

    field[i][0] += (-cphi[i][1] + term * dix);
    field[i][1] += (-cphi[i][2] + term * diy);
    field[i][2] += (-cphi[i][3] + term * diz);
  }
}

// see also subroutine udirect2b / dfield0c in induce.f
void dfield_ewald_real(real* gpu_field, real* gpu_fieldp) {
  real(*field)[3] = reinterpret_cast<real(*)[3]>(gpu_field);
  real(*fieldp)[3] = reinterpret_cast<real(*)[3]>(gpu_fieldp);

  const real off = ewald_switch_off;
  const real off2 = off * off;
  const int maxnlst = mlist_obj_.maxnlst;

  static std::vector<real> pscalebuf;
  static std::vector<real> dscalebuf;
  pscalebuf.resize(n, 1);
  dscalebuf.resize(n, 1);
  real* pscale = pscalebuf.data();
  real* dscale = dscalebuf.data();

  const int pu = ppme_unit;
  const real aewald = pme_obj(pu).aewald;
  const real aesq2 = 2 * aewald * aewald;
  const real aesq2n = (aewald > 0 ? 1 / (sqrtpi * aewald) : 0);

  real bn[4];

  #pragma acc parallel loop independent\
              deviceptr(x,y,z,box,coupl,polargroup,mlst,\
              rpole,thole,pdamp,\
              field,fieldp)\
              firstprivate(pscale[0:n],dscale[0:n])
  for (int i = 0; i < n; ++i) {

    // set exclusion coefficients for connected atoms

    const int n12i = coupl->n12[i];
    const int n13i = coupl->n13[i];
    const int n14i = coupl->n14[i];
    const int n15i = coupl->n15[i];

    const int np11i = polargroup->np11[i];
    const int np12i = polargroup->np12[i];
    const int np13i = polargroup->np13[i];
    const int np14i = polargroup->np14[i];

    #pragma acc loop independent
    for (int j = 0; j < n12i; ++j) {
      int ij = coupl->i12[i][j];
      pscale[ij] = p2scale;
      #pragma acc loop independent
      for (int k = 0; k < np11i; ++k)
        if (ij == polargroup->ip11[i][k])
          pscale[ij] = p2iscale;
    }
    #pragma acc loop independent
    for (int j = 0; j < n13i; ++j) {
      int ij = coupl->i13[i][j];
      pscale[ij] = p3scale;
      #pragma acc loop independent
      for (int k = 0; k < np11i; ++k)
        if (ij == polargroup->ip11[i][k])
          pscale[ij] = p3iscale;
    }
    #pragma acc loop independent
    for (int j = 0; j < n14i; ++j) {
      int ij = coupl->i14[i][j];
      pscale[ij] = p4scale;
      #pragma acc loop independent
      for (int k = 0; k < np11i; ++k)
        if (ij == polargroup->ip11[i][k])
          pscale[ij] = p4iscale;
    }
    #pragma acc loop independent
    for (int j = 0; j < n15i; ++j) {
      int ij = coupl->i15[i][j];
      pscale[ij] = p5scale;
      #pragma acc loop independent
      for (int k = 0; k < np11i; ++k)
        if (ij == polargroup->ip11[i][k])
          pscale[ij] = p5iscale;
    }

    #pragma acc loop independent
    for (int j = 0; j < np11i; ++j)
      dscale[polargroup->ip11[i][j]] = d1scale;
    #pragma acc loop independent
    for (int j = 0; j < np12i; ++j)
      dscale[polargroup->ip12[i][j]] = d2scale;
    #pragma acc loop independent
    for (int j = 0; j < np13i; ++j)
      dscale[polargroup->ip13[i][j]] = d3scale;
    #pragma acc loop independent
    for (int j = 0; j < np14i; ++j)
      dscale[polargroup->ip14[i][j]] = d4scale;

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

    int nmlsti = mlst->nlst[i];
    int base = i * maxnlst;
    #pragma acc loop independent private(bn[0:4])
    for (int kk = 0; kk < nmlsti; ++kk) {
      int k = mlst->lst[base + kk];
      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;

      image(xr, yr, zr, box);
      real r2 = xr * xr + yr * yr + zr * zr;
      if (r2 <= off2) {
        real r = REAL_SQRT(r2);
        real rr1 = REAL_RECIP(r);
        real rr2 = rr1 * rr1;
        real rr3 = rr1 * rr2;
        real rr5 = 3 * rr3 * rr2;
        real rr7 = 5 * rr5 * rr2;
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

        real ralpha = aewald * r;
        bn[0] = REAL_ERFC(ralpha) * rr1;
        real exp2a = REAL_EXP(-REAL_SQ(ralpha));
        real aefac = aesq2n;
        #pragma acc loop seq
        for (int j = 1; j <= 3; ++j) {
          aefac *= aesq2;
          bn[j] = ((j + j - 1) * bn[j - 1] + aefac * exp2a) * rr2;
        }

        // if use_thole
        real scale3 = 1;
        real scale5 = 1;
        real scale7 = 1;
        real damp = pdi * pdamp[k];
        if (damp != 0) {
          real pgamma = REAL_MIN(pti, thole[k]);
          damp = -pgamma * REAL_CUBE(r / damp);
          if (damp > -50) {
            real expdamp = REAL_EXP(damp);
            scale3 = 1 - expdamp;
            scale5 = 1 - expdamp * (1 - damp);
            scale7 = 1 - expdamp * (1 - damp + (real)0.6 * REAL_SQ(damp));
          }
        }
        real scalek, bcn1, bcn2, bcn3;

        scalek = dscale[k];
        bcn1 = bn[1] - (1 - scalek * scale3) * rr3;
        bcn2 = bn[2] - (1 - scalek * scale5) * rr5;
        bcn3 = bn[3] - (1 - scalek * scale7) * rr7;

        real fid1, fid2, fid3;
        real fkd1, fkd2, fkd3;
        fid1 = -xr * (bcn1 * ck - bcn2 * dkr + bcn3 * qkr) - bcn1 * dkx +
            2 * bcn2 * qkx;
        fid2 = -yr * (bcn1 * ck - bcn2 * dkr + bcn3 * qkr) - bcn1 * dky +
            2 * bcn2 * qky;
        fid3 = -zr * (bcn1 * ck - bcn2 * dkr + bcn3 * qkr) - bcn1 * dkz +
            2 * bcn2 * qkz;
        fkd1 = xr * (bcn1 * ci + bcn2 * dir + bcn3 * qir) - bcn1 * dix -
            2 * bcn2 * qix;
        fkd2 = yr * (bcn1 * ci + bcn2 * dir + bcn3 * qir) - bcn1 * diy -
            2 * bcn2 * qiy;
        fkd3 = zr * (bcn1 * ci + bcn2 * dir + bcn3 * qir) - bcn1 * diz -
            2 * bcn2 * qiz;

        scalek = pscale[k];
        bcn1 = bn[1] - (1 - scalek * scale3) * rr3;
        bcn2 = bn[2] - (1 - scalek * scale5) * rr5;
        bcn3 = bn[3] - (1 - scalek * scale7) * rr7;

        real fip1, fip2, fip3;
        real fkp1, fkp2, fkp3;
        fip1 = -xr * (bcn1 * ck - bcn2 * dkr + bcn3 * qkr) - bcn1 * dkx +
            2 * bcn2 * qkx;
        fip2 = -yr * (bcn1 * ck - bcn2 * dkr + bcn3 * qkr) - bcn1 * dky +
            2 * bcn2 * qky;
        fip3 = -zr * (bcn1 * ck - bcn2 * dkr + bcn3 * qkr) - bcn1 * dkz +
            2 * bcn2 * qkz;
        fkp1 = xr * (bcn1 * ci + bcn2 * dir + bcn3 * qir) - bcn1 * dix -
            2 * bcn2 * qix;
        fkp2 = yr * (bcn1 * ci + bcn2 * dir + bcn3 * qir) - bcn1 * diy -
            2 * bcn2 * qiy;
        fkp3 = zr * (bcn1 * ci + bcn2 * dir + bcn3 * qir) - bcn1 * diz -
            2 * bcn2 * qiz;
        // end if use_thole

        #pragma acc atomic update
        field[i][0] += fid1;
        #pragma acc atomic update
        field[i][1] += fid2;
        #pragma acc atomic update
        field[i][2] += fid3;
        #pragma acc atomic update
        field[k][0] += fkd1;
        #pragma acc atomic update
        field[k][1] += fkd2;
        #pragma acc atomic update
        field[k][2] += fkd3;

        #pragma acc atomic update
        fieldp[i][0] += fip1;
        #pragma acc atomic update
        fieldp[i][1] += fip2;
        #pragma acc atomic update
        fieldp[i][2] += fip3;
        #pragma acc atomic update
        fieldp[k][0] += fkp1;
        #pragma acc atomic update
        fieldp[k][1] += fkp2;
        #pragma acc atomic update
        fieldp[k][2] += fkp3;
      }
    } // end for (int kk)

    // reset exclusion coefficients for connected atoms

    #pragma acc loop independent
    for (int j = 0; j < n12i; ++j)
      pscale[coupl->i12[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < n13i; ++j)
      pscale[coupl->i13[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < n14i; ++j)
      pscale[coupl->i14[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < n15i; ++j)
      pscale[coupl->i15[i][j]] = 1;

    #pragma acc loop independent
    for (int j = 0; j < np11i; ++j)
      dscale[polargroup->ip11[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < np12i; ++j)
      dscale[polargroup->ip12[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < np13i; ++j)
      dscale[polargroup->ip13[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < np14i; ++j)
      dscale[polargroup->ip14[i][j]] = 1;
  } // end for (int i)
}

void dfield_ewald(real* gpu_field, real* gpu_fieldp) {
  zero_array(gpu_field, 3 * n);
  zero_array(gpu_fieldp, 3 * n);

  dfield_ewald_recip_self(gpu_field);
  #pragma acc parallel loop independent deviceptr(gpu_field, gpu_fieldp)
  for (int i = 0; i < 3 * n; ++i) {
    gpu_fieldp[i] = gpu_field[i];
  }

  dfield_ewald_real(gpu_field, gpu_fieldp);
}

// see also subroutine umutual1 in induce.f
void ufield_ewald_recip_self(const real* gpu_uind, const real* gpu_uinp,
                             real* gpu_field, real* gpu_fieldp) {
  const real(*uind)[3] = reinterpret_cast<const real(*)[3]>(gpu_uind);
  const real(*uinp)[3] = reinterpret_cast<const real(*)[3]>(gpu_uinp);
  real(*field)[3] = reinterpret_cast<real(*)[3]>(gpu_field);
  real(*fieldp)[3] = reinterpret_cast<real(*)[3]>(gpu_fieldp);

  const int pu = ppme_unit;
  const auto& st = pme_obj(pu);
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

  real a[3][3];

  const real term = REAL_CUBE(aewald) * 4 / 3 / sqrtpi;

  #pragma acc parallel loop independent deviceptr(box,\
              field,fieldp,uind,uinp,\
              fdip_phi1,fdip_phi2)\
              private(a[0:3][0:3])
  for (int i = 0; i < n; ++i) {
    a[0][0] = nfft1 * box->recip[0][0];
    a[1][0] = nfft2 * box->recip[1][0];
    a[2][0] = nfft3 * box->recip[2][0];
    a[0][1] = nfft1 * box->recip[0][1];
    a[1][1] = nfft2 * box->recip[1][1];
    a[2][1] = nfft3 * box->recip[2][1];
    a[0][2] = nfft1 * box->recip[0][2];
    a[1][2] = nfft2 * box->recip[1][2];
    a[2][2] = nfft3 * box->recip[2][2];

    #pragma acc loop independent
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

void ufield_ewald_real(const real* gpu_uind, const real* gpu_uinp,
                       real* gpu_field, real* gpu_fieldp) {
  const real(*uind)[3] = reinterpret_cast<const real(*)[3]>(gpu_uind);
  const real(*uinp)[3] = reinterpret_cast<const real(*)[3]>(gpu_uinp);
  real(*field)[3] = reinterpret_cast<real(*)[3]>(gpu_field);
  real(*fieldp)[3] = reinterpret_cast<real(*)[3]>(gpu_fieldp);

  const real off = ewald_switch_cut;
  const real off2 = off * off;
  const int maxnlst = mlist_obj_.maxnlst;

  static std::vector<real> uscalebuf;
  uscalebuf.resize(n, 1);
  real* uscale = uscalebuf.data();

  const int pu = ppme_unit;
  const real aewald = pme_obj(pu).aewald;
  const real aesq2 = 2 * aewald * aewald;
  const real aesq2n = (aewald > 0 ? 1 / (sqrtpi * aewald) : 0);

  real bn[3];

  #pragma acc parallel loop independent\
              deviceptr(x,y,z,box,polargroup,mlst,\
              thole,pdamp,uind,uinp,field,fieldp)\
              firstprivate(uscale[0:n])
  for (int i = 0; i < n; ++i) {

    // set exclusion coefficients for connected atoms

    const int np11i = polargroup->np11[i];
    const int np12i = polargroup->np12[i];
    const int np13i = polargroup->np13[i];
    const int np14i = polargroup->np14[i];

    #pragma acc loop independent
    for (int j = 0; j < np11i; ++j)
      uscale[polargroup->ip11[i][j]] = u1scale;
    #pragma acc loop independent
    for (int j = 0; j < np12i; ++j)
      uscale[polargroup->ip12[i][j]] = u2scale;
    #pragma acc loop independent
    for (int j = 0; j < np13i; ++j)
      uscale[polargroup->ip13[i][j]] = u3scale;
    #pragma acc loop independent
    for (int j = 0; j < np14i; ++j)
      uscale[polargroup->ip14[i][j]] = u4scale;

    real xi = x[i];
    real yi = y[i];
    real zi = z[i];

    int nmlsti = mlst->nlst[i];
    int base = i * maxnlst;
    #pragma acc loop independent private(bn[0:3])
    for (int kk = 0; kk < nmlsti; ++kk) {
      int k = mlst->lst[base + kk];
      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;
      real pdi = pdamp[i];
      real pti = thole[i];

      image(xr, yr, zr, box);
      real r2 = xr * xr + yr * yr + zr * zr;
      if (r2 <= off2) {
        real r = REAL_SQRT(r2);
        real rr1 = REAL_RECIP(r);
        real rr2 = rr1 * rr1;
        real rr3 = rr1 * rr2;
        real rr5 = 3 * rr3 * rr2;

        real ralpha = aewald * r;
        bn[0] = REAL_ERFC(ralpha) * rr1;
        real exp2a = REAL_EXP(-REAL_SQ(ralpha));
        real aefac = aesq2n;
        #pragma acc loop seq
        for (int j = 1; j <= 2; ++j) {
          aefac *= aesq2;
          bn[j] = ((j + j - 1) * bn[j - 1] + aefac * exp2a) * rr2;
        }

        // if use_thole
        real scale3 = 1;
        real scale5 = 1;
        real damp = pdi * pdamp[k];
        if (damp != 0) {
          real pgamma = REAL_MIN(pti, thole[k]);
          damp = -pgamma * REAL_CUBE(r / damp);
          if (damp > -50) {
            real expdamp = REAL_EXP(damp);
            scale3 = 1 - expdamp;
            scale5 = 1 - expdamp * (1 - damp);
          }
        }
        real scalek, bcn1, bcn2;

        scalek = uscale[k];
        bcn1 = bn[1] - (1 - scalek * scale3) * rr3;
        bcn2 = bn[2] - (1 - scalek * scale5) * rr5;

        real dlocal1 = bcn2 * xr * xr - bcn1;
        real dlocal2 = bcn2 * xr * yr;
        real dlocal3 = bcn2 * xr * zr;
        real dlocal4 = bcn2 * yr * yr - bcn1;
        real dlocal5 = bcn2 * yr * zr;
        real dlocal6 = bcn2 * zr * zr - bcn1;

        real fid1 =
            dlocal1 * uind[k][0] + dlocal2 * uind[k][1] + dlocal3 * uind[k][2];
        real fid2 =
            dlocal2 * uind[k][0] + dlocal4 * uind[k][1] + dlocal5 * uind[k][2];
        real fid3 =
            dlocal3 * uind[k][0] + dlocal5 * uind[k][1] + dlocal6 * uind[k][2];
        real fkd1 =
            dlocal1 * uind[i][0] + dlocal2 * uind[i][1] + dlocal3 * uind[i][2];
        real fkd2 =
            dlocal2 * uind[i][0] + dlocal4 * uind[i][1] + dlocal5 * uind[i][2];
        real fkd3 =
            dlocal3 * uind[i][0] + dlocal5 * uind[i][1] + dlocal6 * uind[i][2];

        real fip1 =
            dlocal1 * uinp[k][0] + dlocal2 * uinp[k][1] + dlocal3 * uinp[k][2];
        real fip2 =
            dlocal2 * uinp[k][0] + dlocal4 * uinp[k][1] + dlocal5 * uinp[k][2];
        real fip3 =
            dlocal3 * uinp[k][0] + dlocal5 * uinp[k][1] + dlocal6 * uinp[k][2];
        real fkp1 =
            dlocal1 * uinp[i][0] + dlocal2 * uinp[i][1] + dlocal3 * uinp[i][2];
        real fkp2 =
            dlocal2 * uinp[i][0] + dlocal4 * uinp[i][1] + dlocal5 * uinp[i][2];
        real fkp3 =
            dlocal3 * uinp[i][0] + dlocal5 * uinp[i][1] + dlocal6 * uinp[i][2];

        #pragma acc atomic update
        field[i][0] += fid1;
        #pragma acc atomic update
        field[i][1] += fid2;
        #pragma acc atomic update
        field[i][2] += fid3;
        #pragma acc atomic update
        field[k][0] += fkd1;
        #pragma acc atomic update
        field[k][1] += fkd2;
        #pragma acc atomic update
        field[k][2] += fkd3;

        #pragma acc atomic update
        fieldp[i][0] += fip1;
        #pragma acc atomic update
        fieldp[i][1] += fip2;
        #pragma acc atomic update
        fieldp[i][2] += fip3;
        #pragma acc atomic update
        fieldp[k][0] += fkp1;
        #pragma acc atomic update
        fieldp[k][1] += fkp2;
        #pragma acc atomic update
        fieldp[k][2] += fkp3;
        // end if use_thole
      }
    } // end for (int kk)

    // reset exclusion coefficients for connected atoms

    #pragma acc loop independent
    for (int j = 0; j < np11i; ++j)
      uscale[polargroup->ip11[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < np12i; ++j)
      uscale[polargroup->ip12[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < np13i; ++j)
      uscale[polargroup->ip13[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < np14i; ++j)
      uscale[polargroup->ip14[i][j]] = 1;
  } // end for (int i)
}

void ufield_ewald(const real* gpu_uind, const real* gpu_uinp, real* gpu_field,
                  real* gpu_fieldp) {
  zero_array(gpu_field, 3 * n);
  zero_array(gpu_fieldp, 3 * n);

  ufield_ewald_recip_self(gpu_uind, gpu_uinp, gpu_field, gpu_fieldp);
  ufield_ewald_real(gpu_uind, gpu_uinp, gpu_field, gpu_fieldp);
}
}
TINKER_NAMESPACE_END
