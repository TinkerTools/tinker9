#include "acc_seq.h"
#include "gpu/e_polar.h"
#include "util_mdstate.h"
#include "util_pme.h"

TINKER_NAMESPACE_BEGIN
template <int USE>
void epolar_real_tmpl(const real (*gpu_uind)[3], const real (*gpu_uinp)[3]) {
  constexpr int do_e = USE & calc::energy;
  constexpr int do_a = USE & calc::analyz;
  constexpr int do_g = USE & calc::grad;
  constexpr int do_v = USE & calc::virial;
  sanity_check<USE>();

  if_constexpr(do_g) {
    zero_array(&ufld[0][0], 3 * n);
    zero_array(&dufld[0][0], 6 * n);
  }

  const real(*uind)[3] = reinterpret_cast<const real(*)[3]>(gpu_uind);
  const real(*uinp)[3] = reinterpret_cast<const real(*)[3]>(gpu_uinp);

  const real off = ewald_switch_off;
  const real off2 = off * off;
  const int maxnlst = mlist_obj_.maxnlst;

  static std::vector<real> pscalebuf;
  static std::vector<real> dscalebuf;
  static std::vector<real> uscalebuf;
  pscalebuf.resize(n, 1);
  dscalebuf.resize(n, 1);
  uscalebuf.resize(n, 1);
  real* pscale = pscalebuf.data();
  real* dscale = dscalebuf.data();
  real* uscale = uscalebuf.data();

  const real f = 0.5 * electric / dielec;

  const real aewald = pme_obj(ppme_unit).aewald;
  real bn[5];

  #pragma acc parallel loop independent\
              deviceptr(x,y,z,box,coupl,polargroup,mlst,\
              rpole,thole,pdamp,uind,uinp,\
              ep,nep,vir_ep,ufld,dufld)\
              firstprivate(pscale[0:n],dscale[0:n],uscale[0:n])
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
          for (int j = 1; j <= 3; ++j) {
            alsq2n *= alsq2;
            bn[j] = ((j + j - 1) * bn[j - 1] + alsq2n * exp2a) * rr2;
          }
        }
        else {
          #pragma acc loop seq
          for (int j = 1; j <= 4; ++j) {
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

        if_constexpr(do_g) {
          real uirp = uixp * xr + uiyp * yr + uizp * zr;
          real ukrp = ukxp * xr + ukyp * yr + ukzp * zr;

          real psc3 = 1 - sc3 * pscale[k];
          real psc5 = 1 - sc5 * pscale[k];
          real psc7 = 1 - sc7 * pscale[k];
          real dsc3 = 1 - sc3 * dscale[k];
          real dsc5 = 1 - sc5 * dscale[k];
          real dsc7 = 1 - sc7 * dscale[k];
          real usc3 = 1 - sc3 * uscale[k];
          real usc5 = 1 - sc5 * uscale[k];
          real psr3 = bn[1] - psc3 * rr3;
          real psr5 = bn[2] - psc5 * rr5;
          real psr7 = bn[3] - psc7 * rr7;
          real dsr3 = bn[1] - dsc3 * rr3;
          real dsr5 = bn[2] - dsc5 * rr5;
          real dsr7 = bn[3] - dsc7 * rr7;
          // variable "usr3" was declared but never referenced
          // real usr3 = bn[1] - usc3 * rr3;
          real usr5 = bn[2] - usc5 * rr5;

          real prc31 = rc31 * pscale[k];
          real prc51 = rc51 * pscale[k];
          real prc71 = rc71 * pscale[k];
          real drc31 = rc31 * dscale[k];
          real drc51 = rc51 * dscale[k];
          real drc71 = rc71 * dscale[k];
          real urc31 = rc31 * uscale[k];
          real urc51 = rc51 * uscale[k];

          real prc32 = rc32 * pscale[k];
          real prc52 = rc52 * pscale[k];
          real prc72 = rc72 * pscale[k];
          real drc32 = rc32 * dscale[k];
          real drc52 = rc52 * dscale[k];
          real drc72 = rc72 * dscale[k];
          real urc32 = rc32 * uscale[k];
          real urc52 = rc52 * uscale[k];

          real prc33 = rc33 * pscale[k];
          real prc53 = rc53 * pscale[k];
          real prc73 = rc73 * pscale[k];
          real drc33 = rc33 * dscale[k];
          real drc53 = rc53 * dscale[k];
          real drc73 = rc73 * dscale[k];
          // variable "urc33" was declared but never referenced
          // real urc33 = rc33 * uscale[k];
          real urc53 = rc53 * uscale[k];

          // get the induced dipole field used for dipole torques

          real tuir, tukr;
          real tix3 = psr3 * ukx + dsr3 * ukxp;
          real tiy3 = psr3 * uky + dsr3 * ukyp;
          real tiz3 = psr3 * ukz + dsr3 * ukzp;
          real tkx3 = psr3 * uix + dsr3 * uixp;
          real tky3 = psr3 * uiy + dsr3 * uiyp;
          real tkz3 = psr3 * uiz + dsr3 * uizp;
          tuir = -psr5 * ukr - dsr5 * ukrp;
          tukr = -psr5 * uir - dsr5 * uirp;

          #pragma acc atomic update
          ufld[i][0] += tix3 + xr * tuir;
          #pragma acc atomic update
          ufld[i][1] += tiy3 + yr * tuir;
          #pragma acc atomic update
          ufld[i][2] += tiz3 + zr * tuir;
          #pragma acc atomic update
          ufld[k][0] += tkx3 + xr * tukr;
          #pragma acc atomic update
          ufld[k][1] += tky3 + yr * tukr;
          #pragma acc atomic update
          ufld[k][2] += tkz3 + zr * tukr;

          // get induced dipole field gradient used for quadrupole torques

          real tix5 = 2 * (psr5 * ukx + dsr5 * ukxp);
          real tiy5 = 2 * (psr5 * uky + dsr5 * ukyp);
          real tiz5 = 2 * (psr5 * ukz + dsr5 * ukzp);
          real tkx5 = 2 * (psr5 * uix + dsr5 * uixp);
          real tky5 = 2 * (psr5 * uiy + dsr5 * uiyp);
          real tkz5 = 2 * (psr5 * uiz + dsr5 * uizp);
          tuir = -psr7 * ukr - dsr7 * ukrp;
          tukr = -psr7 * uir - dsr7 * uirp;

          #pragma acc atomic update
          dufld[i][0] += (xr * tix5 + xr * xr * tuir);
          #pragma acc atomic update
          dufld[i][1] += (xr * tiy5 + yr * tix5 + 2 * xr * yr * tuir);
          #pragma acc atomic update
          dufld[i][2] += (yr * tiy5 + yr * yr * tuir);
          #pragma acc atomic update
          dufld[i][3] += (xr * tiz5 + zr * tix5 + 2 * xr * zr * tuir);
          #pragma acc atomic update
          dufld[i][4] += (yr * tiz5 + zr * tiy5 + 2 * yr * zr * tuir);
          #pragma acc atomic update
          dufld[i][5] += (zr * tiz5 + zr * zr * tuir);

          #pragma acc atomic update
          dufld[k][0] += (-xr * tkx5 - xr * xr * tukr);
          #pragma acc atomic update
          dufld[k][1] += (-xr * tky5 - yr * tkx5 - 2 * xr * yr * tukr);
          #pragma acc atomic update
          dufld[k][2] += (-yr * tky5 - yr * yr * tukr);
          #pragma acc atomic update
          dufld[k][3] += (-xr * tkz5 - zr * tkx5 - 2 * xr * zr * tukr);
          #pragma acc atomic update
          dufld[k][4] += (-yr * tkz5 - zr * tky5 - 2 * yr * zr * tukr);
          #pragma acc atomic update
          dufld[k][5] += (-zr * tkz5 - zr * zr * tukr);

          real term1, term2, term3, term4, term5, term6, term7;
          real depx, depy, depz;
          real frcx, frcy, frcz;

          // get the dEd/dR terms used for direct polarization force

          term1 = bn[2] - dsc3 * rr5;
          term2 = bn[3] - dsc5 * rr7;
          term3 = -dsr3 + term1 * xr * xr - rr3 * xr * drc31;
          term4 = rr3 * drc31 - term1 * xr - dsr5 * xr;
          term5 = term2 * xr * xr - dsr5 - rr5 * xr * drc51;
          term6 = (bn[4] - dsc7 * rr9) * xr * xr - bn[3] - rr7 * xr * drc71;
          term7 =
              rr5 * drc51 - 2 * bn[3] * xr + (dsc5 + 1.5f * dsc7) * rr7 * xr;
          real tixx = ci * term3 + dix * term4 + dir * term5 + 2 * dsr5 * qixx +
              (qiy * yr + qiz * zr) * dsc7 * rr7 + 2 * qix * term7 +
              qir * term6;
          real tkxx = ck * term3 - dkx * term4 - dkr * term5 + 2 * dsr5 * qkxx +
              (qky * yr + qkz * zr) * dsc7 * rr7 + 2 * qkx * term7 +
              qkr * term6;

          term3 = -dsr3 + term1 * yr * yr - rr3 * yr * drc32;
          term4 = rr3 * drc32 - term1 * yr - dsr5 * yr;
          term5 = term2 * yr * yr - dsr5 - rr5 * yr * drc52;
          term6 = (bn[4] - dsc7 * rr9) * yr * yr - bn[3] - rr7 * yr * drc72;
          term7 =
              rr5 * drc52 - 2 * bn[3] * yr + (dsc5 + 1.5f * dsc7) * rr7 * yr;
          real tiyy = ci * term3 + diy * term4 + dir * term5 + 2 * dsr5 * qiyy +
              (qix * xr + qiz * zr) * dsc7 * rr7 + 2 * qiy * term7 +
              qir * term6;
          real tkyy = ck * term3 - dky * term4 - dkr * term5 + 2 * dsr5 * qkyy +
              (qkx * xr + qkz * zr) * dsc7 * rr7 + 2 * qky * term7 +
              qkr * term6;

          term3 = -dsr3 + term1 * zr * zr - rr3 * zr * drc33;
          term4 = rr3 * drc33 - term1 * zr - dsr5 * zr;
          term5 = term2 * zr * zr - dsr5 - rr5 * zr * drc53;
          term6 = (bn[4] - dsc7 * rr9) * zr * zr - bn[3] - rr7 * zr * drc73;
          term7 =
              rr5 * drc53 - 2 * bn[3] * zr + (dsc5 + 1.5f * dsc7) * rr7 * zr;
          real tizz = ci * term3 + diz * term4 + dir * term5 + 2 * dsr5 * qizz +
              (qix * xr + qiy * yr) * dsc7 * rr7 + 2 * qiz * term7 +
              qir * term6;
          real tkzz = ck * term3 - dkz * term4 - dkr * term5 + 2 * dsr5 * qkzz +
              (qkx * xr + qky * yr) * dsc7 * rr7 + 2 * qkz * term7 +
              qkr * term6;

          term3 = term1 * xr * yr - rr3 * yr * drc31;
          term4 = rr3 * drc31 - term1 * xr;
          term5 = term2 * xr * yr - rr5 * yr * drc51;
          term6 = (bn[4] - dsc7 * rr9) * xr * yr - rr7 * yr * drc71;
          term7 = rr5 * drc51 - term2 * xr;
          real tixy = ci * term3 - dsr5 * dix * yr + diy * term4 + dir * term5 +
              2 * dsr5 * qixy - 2 * dsr7 * yr * qix + 2 * qiy * term7 +
              qir * term6;
          real tkxy = ck * term3 + dsr5 * dkx * yr - dky * term4 - dkr * term5 +
              2 * dsr5 * qkxy - 2 * dsr7 * yr * qkx + 2 * qky * term7 +
              qkr * term6;

          term3 = term1 * xr * zr - rr3 * zr * drc31;
          term5 = term2 * xr * zr - rr5 * zr * drc51;
          term6 = (bn[4] - dsc7 * rr9) * xr * zr - rr7 * zr * drc71;
          real tixz = ci * term3 - dsr5 * dix * zr + diz * term4 + dir * term5 +
              2 * dsr5 * qixz - 2 * dsr7 * zr * qix + 2 * qiz * term7 +
              qir * term6;
          real tkxz = ck * term3 + dsr5 * dkx * zr - dkz * term4 - dkr * term5 +
              2 * dsr5 * qkxz - 2 * dsr7 * zr * qkx + 2 * qkz * term7 +
              qkr * term6;

          term3 = term1 * yr * zr - rr3 * zr * drc32;
          term4 = rr3 * drc32 - term1 * yr;
          term5 = term2 * yr * zr - rr5 * zr * drc52;
          term6 = (bn[4] - dsc7 * rr9) * yr * zr - rr7 * zr * drc72;
          term7 = rr5 * drc52 - term2 * yr;
          real tiyz = ci * term3 - dsr5 * diy * zr + diz * term4 + dir * term5 +
              2 * dsr5 * qiyz - 2 * dsr7 * zr * qiy + 2 * qiz * term7 +
              qir * term6;
          real tkyz = ck * term3 + dsr5 * dky * zr - dkz * term4 - dkr * term5 +
              2 * dsr5 * qkyz - 2 * dsr7 * zr * qky + 2 * qkz * term7 +
              qkr * term6;

          depx = tixx * ukxp + tixy * ukyp + tixz * ukzp - tkxx * uixp -
              tkxy * uiyp - tkxz * uizp;
          depy = tixy * ukxp + tiyy * ukyp + tiyz * ukzp - tkxy * uixp -
              tkyy * uiyp - tkyz * uizp;
          depz = tixz * ukxp + tiyz * ukyp + tizz * ukzp - tkxz * uixp -
              tkyz * uiyp - tkzz * uizp;

          frcx = depx;
          frcy = depy;
          frcz = depz;

          // get the dEp/dR terms used for direct polarization force

          term1 = bn[2] - psc3 * rr5;
          term2 = bn[3] - psc5 * rr7;
          term3 = -psr3 + term1 * xr * xr - rr3 * xr * prc31;
          term4 = rr3 * prc31 - term1 * xr - psr5 * xr;
          term5 = term2 * xr * xr - psr5 - rr5 * xr * prc51;
          term6 = (bn[4] - psc7 * rr9) * xr * xr - bn[3] - rr7 * xr * prc71;
          term7 =
              rr5 * prc51 - 2 * bn[3] * xr + (psc5 + 1.5f * psc7) * rr7 * xr;
          tixx = ci * term3 + dix * term4 + dir * term5 + 2 * psr5 * qixx +
              (qiy * yr + qiz * zr) * psc7 * rr7 + 2 * qix * term7 +
              qir * term6;
          tkxx = ck * term3 - dkx * term4 - dkr * term5 + 2 * psr5 * qkxx +
              (qky * yr + qkz * zr) * psc7 * rr7 + 2 * qkx * term7 +
              qkr * term6;

          term3 = -psr3 + term1 * yr * yr - rr3 * yr * prc32;
          term4 = rr3 * prc32 - term1 * yr - psr5 * yr;
          term5 = term2 * yr * yr - psr5 - rr5 * yr * prc52;
          term6 = (bn[4] - psc7 * rr9) * yr * yr - bn[3] - rr7 * yr * prc72;
          term7 =
              rr5 * prc52 - 2 * bn[3] * yr + (psc5 + 1.5f * psc7) * rr7 * yr;
          tiyy = ci * term3 + diy * term4 + dir * term5 + 2 * psr5 * qiyy +
              (qix * xr + qiz * zr) * psc7 * rr7 + 2 * qiy * term7 +
              qir * term6;
          tkyy = ck * term3 - dky * term4 - dkr * term5 + 2 * psr5 * qkyy +
              (qkx * xr + qkz * zr) * psc7 * rr7 + 2 * qky * term7 +
              qkr * term6;

          term3 = -psr3 + term1 * zr * zr - rr3 * zr * prc33;
          term4 = rr3 * prc33 - term1 * zr - psr5 * zr;
          term5 = term2 * zr * zr - psr5 - rr5 * zr * prc53;
          term6 = (bn[4] - psc7 * rr9) * zr * zr - bn[3] - rr7 * zr * prc73;
          term7 =
              rr5 * prc53 - 2 * bn[3] * zr + (psc5 + 1.5f * psc7) * rr7 * zr;
          tizz = ci * term3 + diz * term4 + dir * term5 + 2 * psr5 * qizz +
              (qix * xr + qiy * yr) * psc7 * rr7 + 2 * qiz * term7 +
              qir * term6;
          tkzz = ck * term3 - dkz * term4 - dkr * term5 + 2 * psr5 * qkzz +
              (qkx * xr + qky * yr) * psc7 * rr7 + 2 * qkz * term7 +
              qkr * term6;

          term3 = term1 * xr * yr - rr3 * yr * prc31;
          term4 = rr3 * prc31 - term1 * xr;
          term5 = term2 * xr * yr - rr5 * yr * prc51;
          term6 = (bn[4] - psc7 * rr9) * xr * yr - rr7 * yr * prc71;
          term7 = rr5 * prc51 - term2 * xr;
          tixy = ci * term3 - psr5 * dix * yr + diy * term4 + dir * term5 +
              2 * psr5 * qixy - 2 * psr7 * yr * qix + 2 * qiy * term7 +
              qir * term6;
          tkxy = ck * term3 + psr5 * dkx * yr - dky * term4 - dkr * term5 +
              2 * psr5 * qkxy - 2 * psr7 * yr * qkx + 2 * qky * term7 +
              qkr * term6;

          term3 = term1 * xr * zr - rr3 * zr * prc31;
          term5 = term2 * xr * zr - rr5 * zr * prc51;
          term6 = (bn[4] - psc7 * rr9) * xr * zr - rr7 * zr * prc71;
          tixz = ci * term3 - psr5 * dix * zr + diz * term4 + dir * term5 +
              2 * psr5 * qixz - 2 * psr7 * zr * qix + 2 * qiz * term7 +
              qir * term6;
          tkxz = ck * term3 + psr5 * dkx * zr - dkz * term4 - dkr * term5 +
              2 * psr5 * qkxz - 2 * psr7 * zr * qkx + 2 * qkz * term7 +
              qkr * term6;

          term3 = term1 * yr * zr - rr3 * zr * prc32;
          term4 = rr3 * prc32 - term1 * yr;
          term5 = term2 * yr * zr - rr5 * zr * prc52;
          term6 = (bn[4] - psc7 * rr9) * yr * zr - rr7 * zr * prc72;
          term7 = rr5 * prc52 - term2 * yr;
          tiyz = ci * term3 - psr5 * diy * zr + diz * term4 + dir * term5 +
              2 * psr5 * qiyz - 2 * psr7 * zr * qiy + 2 * qiz * term7 +
              qir * term6;
          tkyz = ck * term3 + psr5 * dky * zr - dkz * term4 - dkr * term5 +
              2 * psr5 * qkyz - 2 * psr7 * zr * qky + 2 * qkz * term7 +
              qkr * term6;

          depx = tixx * ukx + tixy * uky + tixz * ukz - tkxx * uix -
              tkxy * uiy - tkxz * uiz;
          depy = tixy * ukx + tiyy * uky + tiyz * ukz - tkxy * uix -
              tkyy * uiy - tkyz * uiz;
          depz = tixz * ukx + tiyz * uky + tizz * ukz - tkxz * uix -
              tkyz * uiy - tkzz * uiz;
          frcx = frcx + depx;
          frcy = frcy + depy;
          frcz = frcz + depz;

          // get the dtau/dr terms used for mutual polarization force

          term1 = bn[2] - usc3 * rr5;
          term2 = bn[3] - usc5 * rr7;
          term3 = usr5 + term1;
          term4 = rr3 * uscale[k];
          term5 = -xr * term3 + rc31 * term4;
          term6 = -usr5 + xr * xr * term2 - rr5 * xr * urc51;
          tixx = uix * term5 + uir * term6;
          tkxx = ukx * term5 + ukr * term6;

          term5 = -yr * term3 + rc32 * term4;
          term6 = -usr5 + yr * yr * term2 - rr5 * yr * urc52;
          tiyy = uiy * term5 + uir * term6;
          tkyy = uky * term5 + ukr * term6;

          term5 = -zr * term3 + rc33 * term4;
          term6 = -usr5 + zr * zr * term2 - rr5 * zr * urc53;
          tizz = uiz * term5 + uir * term6;
          tkzz = ukz * term5 + ukr * term6;

          term4 = -usr5 * yr;
          term5 = -xr * term1 + rr3 * urc31;
          term6 = xr * yr * term2 - rr5 * yr * urc51;
          tixy = uix * term4 + uiy * term5 + uir * term6;
          tkxy = ukx * term4 + uky * term5 + ukr * term6;

          term4 = -usr5 * zr;
          term6 = xr * zr * term2 - rr5 * zr * urc51;
          tixz = uix * term4 + uiz * term5 + uir * term6;
          tkxz = ukx * term4 + ukz * term5 + ukr * term6;

          term5 = -yr * term1 + rr3 * urc32;
          term6 = yr * zr * term2 - rr5 * zr * urc52;
          tiyz = uiy * term4 + uiz * term5 + uir * term6;
          tkyz = uky * term4 + ukz * term5 + ukr * term6;

          depx = tixx * ukxp + tixy * ukyp + tixz * ukzp + tkxx * uixp +
              tkxy * uiyp + tkxz * uizp;
          depy = tixy * ukxp + tiyy * ukyp + tiyz * ukzp + tkxy * uixp +
              tkyy * uiyp + tkyz * uizp;
          depz = tixz * ukxp + tiyz * ukyp + tizz * ukzp + tkxz * uixp +
              tkyz * uiyp + tkzz * uizp;

          frcx = frcx + depx;
          frcy = frcy + depy;
          frcz = frcz + depz;

          #pragma acc atomic update
          gx[i] -= frcx;
          #pragma acc atomic update
          gy[i] -= frcy;
          #pragma acc atomic update
          gz[i] -= frcz;
          #pragma acc atomic update
          gx[k] += frcx;
          #pragma acc atomic update
          gy[k] += frcy;
          #pragma acc atomic update
          gz[k] += frcz;

          // virial

          if_constexpr(do_v) {
            real vxx = xr * frcx;
            real vxy = 0.5f * (yr * frcx + xr * frcy);
            real vxz = 0.5f * (zr * frcx + xr * frcz);
            real vyy = yr * frcy;
            real vyz = 0.5f * (zr * frcy + yr * frcz);
            real vzz = zr * frcz;

            #pragma acc atomic update
            vir_ep[_xx] += vxx;
            #pragma acc atomic update
            vir_ep[_yx] += vxy;
            #pragma acc atomic update
            vir_ep[_zx] += vxz;
            #pragma acc atomic update
            vir_ep[_xy] += vxy;
            #pragma acc atomic update
            vir_ep[_yy] += vyy;
            #pragma acc atomic update
            vir_ep[_zy] += vyz;
            #pragma acc atomic update
            vir_ep[_xz] += vxz;
            #pragma acc atomic update
            vir_ep[_yz] += vyz;
            #pragma acc atomic update
            vir_ep[_zz] += vzz;
          }
        }
        // end if use_thole
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

  // torque

  if_constexpr(do_g) {
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
void epolar_recip_self_tmpl(const real (*gpu_uind)[3],
                            const real (*gpu_uinp)[3]) {
  constexpr int do_e = USE & calc::energy;
  constexpr int do_a = USE & calc::analyz;
  // constexpr int do_g = USE & calc::grad;
  constexpr int do_v = USE & calc::virial;
  sanity_check<USE>();

  const int pu = ppme_unit;
  const auto& st = pme_obj(pu);
  const int nfft1 = st.nfft1;
  const int nfft2 = st.nfft2;
  const int nfft3 = st.nfft3;
  const real aewald = st.aewald;

  const real f = electric / dielec;

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
  #pragma acc parallel loop independent\
              deviceptr(ep,nep,trqx,trqy,trqz,\
              rpole,cmp,gpu_uind,gpu_uinp,cphidp)
  for (int i = 0; i < n; ++i) {
    real tep1 = cmp[i][3] * cphidp[i][2] - cmp[i][2] * cphidp[i][3] +
        2 * (cmp[i][6] - cmp[i][5]) * cphidp[i][9] + cmp[i][8] * cphidp[i][7] +
        cmp[i][9] * cphidp[i][5] - cmp[i][7] * cphidp[i][8] -
        cmp[i][9] * cphidp[i][6];
    real tep2 = cmp[i][1] * cphidp[i][3] - cmp[i][3] * cphidp[i][1] +
        2 * (cmp[i][4] - cmp[i][6]) * cphidp[i][8] + cmp[i][7] * cphidp[i][9] +
        cmp[i][8] * cphidp[i][6] - cmp[i][8] * cphidp[i][4] -
        cmp[i][9] * cphidp[i][7];
    real tep3 = cmp[i][2] * cphidp[i][1] - cmp[i][1] * cphidp[i][2] +
        2 * (cmp[i][5] - cmp[i][4]) * cphidp[i][7] + cmp[i][7] * cphidp[i][4] +
        cmp[i][9] * cphidp[i][8] - cmp[i][7] * cphidp[i][5] -
        cmp[i][8] * cphidp[i][9];

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

  // recip virial

  if_constexpr(do_v) {

    #pragma acc parallel loop independent deviceptr(vir_ep,vir_m)
    for (int i = 0; i < 9; ++i) {
      vir_ep[i] -= vir_m[i];
    }

    scale_data(&cphi[0][0], f, 10 * n);
    scale_data(&fphid[0][0], f, 10 * n);
    scale_data(&fphip[0][0], f, 10 * n);

    real cphid[4], cphip[4];
    real ftc[3][3];
    #pragma acc parallel loop independent deviceptr(vir_ep,box,cmp,\
                gpu_uind,gpu_uinp,fphid,fphip,cphi,cphidp)\
                private(cphid[0:4],cphip[0:4],ftc[0:3][0:3])
    for (int i = 0; i < n; ++i) {

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
      vxy = vxy - 0.5f * (cphidp[i][1] * cmp[i][2] + cphidp[i][2] * cmp[i][1]) -
          0.25f *
              ((gpu_uind[i][1] + gpu_uinp[i][1]) * cphi[i][1] +
               (gpu_uind[i][0] + gpu_uinp[i][0]) * cphi[i][2]);
      vxz = vxz - 0.5f * (cphidp[i][1] * cmp[i][3] + cphidp[i][3] * cmp[i][1]) -
          0.25f *
              ((gpu_uind[i][2] + gpu_uinp[i][2]) * cphi[i][1] +
               (gpu_uind[i][0] + gpu_uinp[i][0]) * cphi[i][3]);
      vyy = vyy - cmp[i][2] * cphidp[i][2] -
          0.5f * ((gpu_uind[i][1] + gpu_uinp[i][1]) * cphi[i][2]);
      vyz = vyz - 0.5f * (cphidp[i][2] * cmp[i][3] + cphidp[i][3] * cmp[i][2]) -
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
      vxx =
          vxx - 0.5f * (cphid[1] * gpu_uinp[i][0] + cphip[1] * gpu_uind[i][0]);
      vxy = vxy -
          0.25f *
              (cphid[1] * gpu_uinp[i][1] + cphip[1] * gpu_uind[i][1] +
               cphid[2] * gpu_uinp[i][0] + cphip[2] * gpu_uind[i][0]);
      vxz = vxz -
          0.25f *
              (cphid[1] * gpu_uinp[i][2] + cphip[1] * gpu_uind[i][2] +
               cphid[3] * gpu_uinp[i][0] + cphip[3] * gpu_uind[i][0]);
      vyy =
          vyy - 0.5f * (cphid[2] * gpu_uinp[i][1] + cphip[2] * gpu_uind[i][1]);
      vyz = vyz -
          0.25f *
              (cphid[2] * gpu_uinp[i][2] + cphip[2] * gpu_uind[i][2] +
               cphid[3] * gpu_uinp[i][1] + cphip[3] * gpu_uind[i][1]);
      vzz =
          vzz - 0.5f * (cphid[3] * gpu_uinp[i][2] + cphip[3] * gpu_uind[i][2]);
      // end if

      #pragma acc atomic update
      vir_ep[_xx] += vxx;
      #pragma acc atomic update
      vir_ep[_yx] += vxy;
      #pragma acc atomic update
      vir_ep[_zx] += vxz;
      #pragma acc atomic update
      vir_ep[_xy] += vxy;
      #pragma acc atomic update
      vir_ep[_yy] += vyy;
      #pragma acc atomic update
      vir_ep[_zy] += vyz;
      #pragma acc atomic update
      vir_ep[_xz] += vxz;
      #pragma acc atomic update
      vir_ep[_yz] += vyz;
      #pragma acc atomic update
      vir_ep[_zz] += vzz;
    }

    // qgrip: pvu_qgrid
    const int pvu = pvpme_unit;
    #pragma acc parallel loop independent deviceptr(cmp,gpu_uinp)
    for (int i = 0; i < n; ++i) {
      cmp[i][1] += gpu_uinp[i][0];
      cmp[i][2] += gpu_uinp[i][1];
      cmp[i][3] += gpu_uinp[i][2];
    }
    cmp_to_fmp(pvu, cmp, fmp);
    grid_mpole(pvu, fmp);
    fftfront(pvu);

    // qgrid: pu_qgrid
    #pragma acc parallel loop independent deviceptr(cmp,gpu_uind,gpu_uinp)
    for (int i = 0; i < n; ++i) {
      cmp[i][1] += (gpu_uind[i][0] - gpu_uinp[i][0]);
      cmp[i][2] += (gpu_uind[i][1] - gpu_uinp[i][1]);
      cmp[i][3] += (gpu_uind[i][2] - gpu_uinp[i][2]);
    }
    cmp_to_fmp(pu, cmp, fmp);
    grid_mpole(pu, fmp);
    fftfront(pu);

    const auto* d = pme_deviceptr(pu);
    const auto* p = pme_deviceptr(pvu);
    const int nff = nfft1 * nfft2;
    const int ntot = nfft1 * nfft2 * nfft3;
    real pterm = REAL_SQ(pi / aewald);

    #pragma acc parallel loop independent deviceptr(box,d,p,vir_ep)
    for (int i = 1; i < ntot; ++i) {
      const real volterm = pi * box->volbox;

      int k3 = i / nff;
      int j = i - k3 * nff;
      int k2 = j / nfft1;
      int k1 = j - k2 * nfft1;

      int r1 = (k1 < (nfft1 + 1) / 2) ? k1 : (k1 - nfft1);
      int r2 = (k2 < (nfft2 + 1) / 2) ? k2 : (k2 - nfft2);
      int r3 = (k3 < (nfft3 + 1) / 2) ? k3 : (k3 - nfft3);

      real h1 =
          box->recip[0][0] * r1 + box->recip[1][0] * r2 + box->recip[2][0] * r3;
      real h2 =
          box->recip[0][1] * r1 + box->recip[1][1] * r2 + box->recip[2][1] * r3;
      real h3 =
          box->recip[0][2] * r1 + box->recip[1][2] * r2 + box->recip[2][2] * r3;
      real hsq = h1 * h1 + h2 * h2 + h3 * h3;
      real term = -pterm * hsq;
      real expterm = 0;
      if (term > -50) {
        // TODO: if .not. use_bounds; if octahedron; 2/hsq
        real denom =
            volterm * hsq * d->bsmod1[k1] * d->bsmod1[k2] * d->bsmod1[k3];
        expterm = REAL_EXP(term) / denom;

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

        #pragma acc atomic update
        vir_ep[_xx] += vxx;
        #pragma acc atomic update
        vir_ep[_xy] += vxy;
        #pragma acc atomic update
        vir_ep[_xz] += vxz;
        #pragma acc atomic update
        vir_ep[_yx] += vxy;
        #pragma acc atomic update
        vir_ep[_yy] += vyy;
        #pragma acc atomic update
        vir_ep[_yz] += vyz;
        #pragma acc atomic update
        vir_ep[_zx] += vxz;
        #pragma acc atomic update
        vir_ep[_zy] += vyz;
        #pragma acc atomic update
        vir_ep[_zz] += vzz;
      }
    }
  }
}

template <int USE>
void epolar_ewald_tmpl(const real (*gpu_uind)[3], const real (*gpu_uinp)[3]) {
  constexpr int do_e = USE & calc::energy;
  constexpr int do_a = USE & calc::analyz;
  constexpr int do_g = USE & calc::grad;
  sanity_check<USE>();

  if_constexpr(do_e && !do_a) epolar0_dotprod(gpu_uind, udirp);
  static_assert(do_g || do_a,
                "Do not use this template for the energy-only version.");

  epolar_real_tmpl<USE>(gpu_uind, gpu_uinp);

  epolar_recip_self_tmpl<USE>(gpu_uind, gpu_uinp);
}

void epolar_ewald(int vers) {
  if (vers == calc::v0) {
    induce(&uind[0][0], &uinp[0][0]);
    epolar0_dotprod(uind, udirp);
  } else if (vers == calc::v1) {
    induce(&uind[0][0], &uinp[0][0]);
    epolar_ewald_tmpl<calc::v1>(uind, uinp);
  } else if (vers == calc::v3) {
    induce(&uind[0][0], &uinp[0][0]);
    epolar_ewald_tmpl<calc::v3>(uind, uinp);
  } else if (vers == calc::v4) {
    induce(&uind[0][0], &uinp[0][0]);
    epolar_ewald_tmpl<calc::v4>(uind, uinp);
  } else if (vers == calc::v5) {
    induce(&uind[0][0], &uinp[0][0]);
    epolar_ewald_tmpl<calc::v5>(uind, uinp);
  } else if (vers == calc::v6) {
    induce(&uind[0][0], &uinp[0][0]);
    epolar_ewald_tmpl<calc::v6>(uind, uinp);
  }
}
TINKER_NAMESPACE_END
