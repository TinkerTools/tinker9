#include "acc_add.h"
#include "acc_image.h"

#include "couple.h"
#include "e_polar.h"
#include "md.h"
#include "nblist.h"

TINKER_NAMESPACE_BEGIN
template <int USE>
void epolar_coulomb_tmpl(const real (*gpu_uind)[3], const real (*gpu_uinp)[3]) {
  constexpr int do_e = USE & calc::energy;
  constexpr int do_a = USE & calc::analyz;
  constexpr int do_g = USE & calc::grad;
  constexpr int do_v = USE & calc::virial;
  static_assert(do_v ? do_g : true, "");
  static_assert(do_a ? do_e : true, "");

  if_constexpr(do_g) { device_array::zero(n, ufld, dufld); }

  if_constexpr(do_e && !do_a) epolar0_dotprod(gpu_uind, udirp);
  static_assert(do_g || do_a,
                "Do not use this template for the energy-only version.");

  const real(*uind)[3] = reinterpret_cast<const real(*)[3]>(gpu_uind);
  const real(*uinp)[3] = reinterpret_cast<const real(*)[3]>(gpu_uinp);

  const real off = mpole_switch_off;
  const real off2 = off * off;
  const int maxnlst = mlist_unit->maxnlst;
  const auto* mlst = mlist_unit.deviceptr();

  const auto* coupl = couple_unit.deviceptr();
  const auto* polargroup = polargroup_unit.deviceptr();

  auto* nep = ep_handle.ne()->buffer();
  auto* ep = ep_handle.e()->buffer();
  auto* vir_ep = ep_handle.vir()->buffer();
  auto bufsize = ep_handle.buffer_size();

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

  #pragma acc parallel num_gangs(bufsize)\
              deviceptr(x,y,z,box,coupl,polargroup,mlst,\
              rpole,thole,pdamp,uind,uinp,\
              ep,nep,vir_ep,ufld,dufld)\
              firstprivate(pscale[0:n],dscale[0:n],uscale[0:n])
  #pragma acc loop gang independent
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
    #pragma acc loop independent
    for (int kk = 0; kk < nmlsti; ++kk) {
      int offset = kk & (bufsize - 1);
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
        real uir = uix * xr + uiy * yr + uiz * zr;
        real ukr = ukx * xr + uky * yr + ukz * zr;

        real r = REAL_SQRT(r2);
        real rr1 = REAL_RECIP(r);
        real rr2 = rr1 * rr1;
        rr1 *= f;
        real rr3 = rr1 * rr2;
        real rr5 = 3 * rr3 * rr2;
        real rr7 = 5 * rr5 * rr2;
        MAYBE_UNUSED real rr9;
        if_constexpr(do_g) { rr9 = 7 * rr7 * rr2; }

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
              real temp3 = -damp * expdamp * rr5;
              real temp5 = -3 * damp * rr2;
              real temp7 = -(1 + 3 * damp) * rr2;
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
        real sr3 = rr3 * sc3;
        real sr5 = rr5 * sc5;
        real sr7 = rr7 * sc7;

        if_constexpr(do_e && do_a) {
          real diu = dix * ukx + diy * uky + diz * ukz;
          real qiu = qix * ukx + qiy * uky + qiz * ukz;
          real dku = dkx * uix + dky * uiy + dkz * uiz;
          real qku = qkx * uix + qky * uiy + qkz * uiz;
          real term1 = ck * uir - ci * ukr + diu + dku;
          real term2 = 2 * (qiu - qku) - uir * dkr - dir * ukr;
          real term3 = uir * qkr - ukr * qir;
          real e = pscale[k] * (term1 * sr3 + term2 * sr5 + term3 * sr7);
          if (e != 0) {
            atomic_add_value(e, ep, offset);
            atomic_add_value(1, nep, offset);
          }
        }

        if_constexpr(do_g) {
          real uirp = uixp * xr + uiyp * yr + uizp * zr;
          real ukrp = ukxp * xr + ukyp * yr + ukzp * zr;

          real dsr3, dsr5, dsr7, psr3, psr5, psr7;
          dsr3 = sr3 * dscale[k];
          dsr5 = sr5 * dscale[k];
          dsr7 = sr7 * dscale[k];
          psr3 = sr3 * pscale[k];
          psr5 = sr5 * pscale[k];
          psr7 = sr7 * pscale[k];

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

          // get the field gradient for direct polarization force

          real term1, term2, term3, term4, term5, term6, term7, term8;

          term1 = sc3 * (rr3 - rr5 * xr * xr) + rc31 * xr;
          term2 = (sc3 + sc5) * rr5 * xr - rc31;
          term3 = sc5 * (rr7 * xr * xr - rr5) - rc51 * xr;
          term4 = 2 * sc5 * rr5;
          term5 = 2 * (sc5 * rr7 * xr - rc51 + 1.5f * sc7 * rr7 * xr);
          term6 = xr * (sc7 * rr9 * xr - rc71);
          real tixx = ci * term1 + dix * term2 - dir * term3 - qixx * term4 +
              qix * term5 - qir * term6 + (qiy * yr + qiz * zr) * sc7 * rr7;
          real tkxx = ck * term1 - dkx * term2 + dkr * term3 - qkxx * term4 +
              qkx * term5 - qkr * term6 + (qky * yr + qkz * zr) * sc7 * rr7;

          term1 = sc3 * (rr3 - rr5 * yr * yr) + rc32 * yr;
          term2 = (sc3 + sc5) * rr5 * yr - rc32;
          term3 = sc5 * (rr7 * yr * yr - rr5) - rc52 * yr;
          term4 = 2 * sc5 * rr5;
          term5 = 2 * (sc5 * rr7 * yr - rc52 + 1.5f * sc7 * rr7 * yr);
          term6 = yr * (sc7 * rr9 * yr - rc72);
          real tiyy = ci * term1 + diy * term2 - dir * term3 - qiyy * term4 +
              qiy * term5 - qir * term6 + (qix * xr + qiz * zr) * sc7 * rr7;
          real tkyy = ck * term1 - dky * term2 + dkr * term3 - qkyy * term4 +
              qky * term5 - qkr * term6 + (qkx * xr + qkz * zr) * sc7 * rr7;

          term1 = sc3 * (rr3 - rr5 * zr * zr) + rc33 * zr;
          term2 = (sc3 + sc5) * rr5 * zr - rc33;
          term3 = sc5 * (rr7 * zr * zr - rr5) - rc53 * zr;
          term4 = 2 * sc5 * rr5;
          term5 = 2 * (sc5 * rr7 * zr - rc53 + 1.5f * sc7 * rr7 * zr);
          term6 = zr * (sc7 * rr9 * zr - rc73);
          real tizz = ci * term1 + diz * term2 - dir * term3 - qizz * term4 +
              qiz * term5 - qir * term6 + (qix * xr + qiy * yr) * sc7 * rr7;
          real tkzz = ck * term1 - dkz * term2 + dkr * term3 - qkzz * term4 +
              qkz * term5 - qkr * term6 + (qkx * xr + qky * yr) * sc7 * rr7;

          term2 = sc3 * rr5 * xr - rc31;
          term1 = yr * term2;
          term3 = sc5 * rr5 * yr;
          term4 = yr * (sc5 * rr7 * xr - rc51);
          term5 = 2 * sc5 * rr5;
          term6 = 2 * (sc5 * rr7 * xr - rc51);
          term7 = 2 * sc7 * rr7 * yr;
          term8 = yr * (sc7 * rr9 * xr - rc71);
          real tixy = -ci * term1 + diy * term2 + dix * term3 - dir * term4 -
              qixy * term5 + qiy * term6 + qix * term7 - qir * term8;
          real tkxy = -ck * term1 - dky * term2 - dkx * term3 + dkr * term4 -
              qkxy * term5 + qky * term6 + qkx * term7 - qkr * term8;

          term2 = sc3 * rr5 * xr - rc31;
          term1 = zr * term2;
          term3 = sc5 * rr5 * zr;
          term4 = zr * (sc5 * rr7 * xr - rc51);
          term5 = 2 * sc5 * rr5;
          term6 = 2 * (sc5 * rr7 * xr - rc51);
          term7 = 2 * sc7 * rr7 * zr;
          term8 = zr * (sc7 * rr9 * xr - rc71);
          real tixz = -ci * term1 + diz * term2 + dix * term3 - dir * term4 -
              qixz * term5 + qiz * term6 + qix * term7 - qir * term8;
          real tkxz = -ck * term1 - dkz * term2 - dkx * term3 + dkr * term4 -
              qkxz * term5 + qkz * term6 + qkx * term7 - qkr * term8;

          term2 = sc3 * rr5 * yr - rc32;
          term1 = zr * term2;
          term3 = sc5 * rr5 * zr;
          term4 = zr * (sc5 * rr7 * yr - rc52);
          term5 = 2 * sc5 * rr5;
          term6 = 2 * (sc5 * rr7 * yr - rc52);
          term7 = 2 * sc7 * rr7 * zr;
          term8 = zr * (sc7 * rr9 * yr - rc72);
          real tiyz = -ci * term1 + diz * term2 + diy * term3 - dir * term4 -
              qiyz * term5 + qiz * term6 + qiy * term7 - qir * term8;
          real tkyz = -ck * term1 - dkz * term2 - dky * term3 + dkr * term4 -
              qkyz * term5 + qkz * term6 + qky * term7 - qkr * term8;

          real depx, depy, depz;
          real frcx, frcy, frcz;

          // get the dEd/dR terms for Thole direct polarization force

          depx = tixx * ukxp + tixy * ukyp + tixz * ukzp - tkxx * uixp -
              tkxy * uiyp - tkxz * uizp;
          depy = tixy * ukxp + tiyy * ukyp + tiyz * ukzp - tkxy * uixp -
              tkyy * uiyp - tkyz * uizp;
          depz = tixz * ukxp + tiyz * ukyp + tizz * ukzp - tkxz * uixp -
              tkyz * uiyp - tkzz * uizp;
          frcx = dscale[k] * depx;
          frcy = dscale[k] * depy;
          frcz = dscale[k] * depz;

          // get the dEp/dR terms for Thole direct polarization force

          depx = tixx * ukx + tixy * uky + tixz * ukz - tkxx * uix -
              tkxy * uiy - tkxz * uiz;
          depy = tixy * ukx + tiyy * uky + tiyz * ukz - tkxy * uix -
              tkyy * uiy - tkyz * uiz;
          depz = tixz * ukx + tiyz * uky + tizz * ukz - tkxz * uix -
              tkyz * uiy - tkzz * uiz;
          frcx += pscale[k] * depx;
          frcy += pscale[k] * depy;
          frcz += pscale[k] * depz;

          // get the dtau/dr terms used for mutual polarization force

          term1 = (sc3 + sc5) * rr5;
          term2 = term1 * xr - rc31;
          term3 = sc5 * (rr5 - rr7 * xr * xr) + rc51 * xr;
          tixx = uix * term2 + uir * term3;
          tkxx = ukx * term2 + ukr * term3;

          term2 = term1 * yr - rc32;
          term3 = sc5 * (rr5 - rr7 * yr * yr) + rc52 * yr;
          tiyy = uiy * term2 + uir * term3;
          tkyy = uky * term2 + ukr * term3;

          term2 = term1 * zr - rc33;
          term3 = sc5 * (rr5 - rr7 * zr * zr) + rc53 * zr;
          tizz = uiz * term2 + uir * term3;
          tkzz = ukz * term2 + ukr * term3;

          term1 = sc5 * rr5 * yr;
          term2 = sc3 * rr5 * xr - rc31;
          term3 = yr * (sc5 * rr7 * xr - rc51);
          tixy = uix * term1 + uiy * term2 - uir * term3;
          tkxy = ukx * term1 + uky * term2 - ukr * term3;
          term1 = sc5 * rr5 * zr;
          term3 = zr * (sc5 * rr7 * xr - rc51);
          tixz = uix * term1 + uiz * term2 - uir * term3;
          tkxz = ukx * term1 + ukz * term2 - ukr * term3;
          term2 = sc3 * rr5 * yr - rc32;
          term3 = zr * (sc5 * rr7 * yr - rc52);
          tiyz = uiy * term1 + uiz * term2 - uir * term3;
          tkyz = uky * term1 + ukz * term2 - ukr * term3;

          depx = tixx * ukxp + tixy * ukyp + tixz * ukzp + tkxx * uixp +
              tkxy * uiyp + tkxz * uizp;
          depy = tixy * ukxp + tiyy * ukyp + tiyz * ukzp + tkxy * uixp +
              tkyy * uiyp + tkyz * uizp;
          depz = tixz * ukxp + tiyz * ukyp + tizz * ukzp + tkxz * uixp +
              tkyz * uiyp + tkzz * uizp;
          frcx += uscale[k] * depx;
          frcy += uscale[k] * depy;
          frcz += uscale[k] * depz;

          #pragma acc atomic update
          gx[i] += frcx;
          #pragma acc atomic update
          gy[i] += frcy;
          #pragma acc atomic update
          gz[i] += frcz;
          #pragma acc atomic update
          gx[k] -= frcx;
          #pragma acc atomic update
          gy[k] -= frcy;
          #pragma acc atomic update
          gz[k] -= frcz;

          // virial

          if_constexpr(do_v) {
            real vxx = -xr * frcx;
            real vxy = -0.5f * (yr * frcx + xr * frcy);
            real vxz = -0.5f * (zr * frcx + xr * frcz);
            real vyy = -yr * frcy;
            real vyz = -0.5f * (zr * frcy + yr * frcz);
            real vzz = -zr * frcz;

            int offv = offset * 16;
            atomic_add_value(vxx, vir_ep, offv + 0);
            atomic_add_value(vxy, vir_ep, offv + 1);
            atomic_add_value(vxz, vir_ep, offv + 2);
            atomic_add_value(vxy, vir_ep, offv + 3);
            atomic_add_value(vyy, vir_ep, offv + 4);
            atomic_add_value(vyz, vir_ep, offv + 5);
            atomic_add_value(vxz, vir_ep, offv + 6);
            atomic_add_value(vyz, vir_ep, offv + 7);
            atomic_add_value(vzz, vir_ep, offv + 8);
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

void epolar_coulomb(int vers) {
  if (vers == calc::v0) {
    induce(&uind[0][0], &uinp[0][0]);
    epolar0_dotprod(uind, udirp);
  } else if (vers == calc::v1) {
    induce(&uind[0][0], &uinp[0][0]);
    epolar_coulomb_tmpl<calc::v1>(uind, uinp);
  } else if (vers == calc::v3) {
    induce(&uind[0][0], &uinp[0][0]);
    epolar_coulomb_tmpl<calc::v3>(uind, uinp);
  } else if (vers == calc::v4) {
    induce(&uind[0][0], &uinp[0][0]);
    epolar_coulomb_tmpl<calc::v4>(uind, uinp);
  } else if (vers == calc::v5) {
    induce(&uind[0][0], &uinp[0][0]);
    epolar_coulomb_tmpl<calc::v5>(uind, uinp);
  } else if (vers == calc::v6) {
    induce(&uind[0][0], &uinp[0][0]);
    epolar_coulomb_tmpl<calc::v6>(uind, uinp);
  }
}
TINKER_NAMESPACE_END
