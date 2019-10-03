#ifndef TINKER_ACC_POLAR_PAIR_H_
#define TINKER_ACC_POLAR_PAIR_H_

#include "acc_damp.h"
#include "acc_mathfunc.h"
#include "e_mpole.h"
#include "e_polar.h"
#include "md.h"

TINKER_NAMESPACE_BEGIN
struct PolarPairGrad {
  real frcx, frcy, frcz;
  real ufldi[3], ufldk[3];
  real dufldi[6], dufldk[6];
};

#pragma acc routine seq
template <int USE, elec_t ETYP>
inline void epolar_pair_acc(                                    //
    real r2, real xr, real yr, real zr,                         //
    real f, real dscale, real pscale, real uscale, real aewald, //
    //
    real ci, real dix, real diy, real diz,                            //
    real qixx, real qixy, real qixz, real qiyy, real qiyz, real qizz, //
    //
    real uix, real uiy, real uiz,                                           //
    MAYBE_UNUSED real uixp, MAYBE_UNUSED real uiyp, MAYBE_UNUSED real uizp, //
    real pdi, real pti,                                                     //
    //
    real ck, real dkx, real dky, real dkz,                            //
    real qkxx, real qkxy, real qkxz, real qkyy, real qkyz, real qkzz, //
    //
    real ukx, real uky, real ukz,                                           //
    MAYBE_UNUSED real ukxp, MAYBE_UNUSED real ukyp, MAYBE_UNUSED real ukzp, //
    real pdk, real ptk,                                                     //
    real& __restrict__ e, PolarPairGrad& __restrict__ pgrad) {
  constexpr int do_e = USE & calc::energy;
  constexpr int do_g = USE & calc::grad;
  constexpr int do_a = USE & calc::analyz;

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
  real invr1 = REAL_RECIP(r);
  real rr2 = invr1 * invr1;

  MAYBE_UNUSED real bn[5];
  if_constexpr(ETYP == elec_t::ewald) {
    if_constexpr(!do_g) damp_ewald(bn, 4, r, invr1, rr2, aewald);
    else damp_ewald(bn, 5, r, invr1, rr2, aewald);
  }
  real rr1 = invr1;
  real rr3 = rr1 * rr2;
  real rr5 = 3 * rr3 * rr2;
  real rr7 = 5 * rr5 * rr2;
  MAYBE_UNUSED real rr9;
  if_constexpr(do_g) rr9 = 7 * rr7 * rr2;

  // if use_thole
  real sc3, sc5, sc7;
  MAYBE_UNUSED real rc31, rc32, rc33, rc51, rc52, rc53, rc71, rc72, rc73;
  if_constexpr(!do_g)                    //
      damp_thole3(r, pdi, pti, pdk, ptk, //
                  sc3, sc5, sc7);
  else                        //
      damp_thole3g(           //
          r, rr2, rr5,        //
          xr, yr, zr,         //
          pdi, pti, pdk, ptk, //
          sc3, sc5, sc7,      //
          rc31, rc32, rc33,   //
          rc51, rc52, rc53,   //
          rc71, rc72, rc73);
  // end if use_thole

  if_constexpr(do_e && do_a) {
    real diu = dix * ukx + diy * uky + diz * ukz;
    real qiu = qix * ukx + qiy * uky + qiz * ukz;
    real dku = dkx * uix + dky * uiy + dkz * uiz;
    real qku = qkx * uix + qky * uiy + qkz * uiz;
    real term1 = ck * uir - ci * ukr + diu + dku;
    real term2 = 2 * (qiu - qku) - uir * dkr - dir * ukr;
    real term3 = uir * qkr - ukr * qir;
    real sr3 = pscale * rr3 * sc3;
    real sr5 = pscale * rr5 * sc5;
    real sr7 = pscale * rr7 * sc7;
    if_constexpr(ETYP == elec_t::ewald) {
      sr3 += (bn[1] - rr3);
      sr5 += (bn[2] - rr5);
      sr7 += (bn[3] - rr7);
    }
    e = f * (term1 * sr3 + term2 * sr5 + term3 * sr7);
  }

  if_constexpr(do_g) {
    real uirp = uixp * xr + uiyp * yr + uizp * zr;
    real ukrp = ukxp * xr + ukyp * yr + ukzp * zr;

    real dsr3, dsr5, dsr7, psr3, psr5, psr7;
    if_constexpr(ETYP == elec_t::ewald) {
      dsr3 = bn[1] - (1 - sc3) * rr3;
      dsr5 = bn[2] - (1 - sc5) * rr5;
      dsr7 = bn[3] - (1 - sc7) * rr7;
    }
    else if_constexpr(ETYP == elec_t::coulomb) {
      dsr3 = rr3 * sc3 * dscale;
      dsr5 = rr5 * sc5 * dscale;
      dsr7 = rr7 * sc7 * dscale;
      psr3 = rr3 * sc3 * pscale;
      psr5 = rr5 * sc5 * pscale;
      psr7 = rr7 * sc7 * pscale;
    }

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

    pgrad.ufldi[0] = f * (tix3 + xr * tuir);
    pgrad.ufldi[1] = f * (tiy3 + yr * tuir);
    pgrad.ufldi[2] = f * (tiz3 + zr * tuir);
    pgrad.ufldk[0] = f * (tkx3 + xr * tukr);
    pgrad.ufldk[1] = f * (tky3 + yr * tukr);
    pgrad.ufldk[2] = f * (tkz3 + zr * tukr);

    // get induced dipole field gradient used for quadrupole torques

    real tix5 = 2 * (psr5 * ukx + dsr5 * ukxp);
    real tiy5 = 2 * (psr5 * uky + dsr5 * ukyp);
    real tiz5 = 2 * (psr5 * ukz + dsr5 * ukzp);
    real tkx5 = 2 * (psr5 * uix + dsr5 * uixp);
    real tky5 = 2 * (psr5 * uiy + dsr5 * uiyp);
    real tkz5 = 2 * (psr5 * uiz + dsr5 * uizp);
    tuir = -psr7 * ukr - dsr7 * ukrp;
    tukr = -psr7 * uir - dsr7 * uirp;

    pgrad.dufldi[0] = f * (xr * tix5 + xr * xr * tuir);
    pgrad.dufldi[1] = f * (xr * tiy5 + yr * tix5 + 2 * xr * yr * tuir);
    pgrad.dufldi[2] = f * (yr * tiy5 + yr * yr * tuir);
    pgrad.dufldi[3] = f * (xr * tiz5 + zr * tix5 + 2 * xr * zr * tuir);
    pgrad.dufldi[4] = f * (yr * tiz5 + zr * tiy5 + 2 * yr * zr * tuir);
    pgrad.dufldi[5] = f * (zr * tiz5 + zr * zr * tuir);
    pgrad.dufldk[0] = f * (-xr * tkx5 - xr * xr * tukr);
    pgrad.dufldk[1] = f * (-xr * tky5 - yr * tkx5 - 2 * xr * yr * tukr);
    pgrad.dufldk[2] = f * (-yr * tky5 - yr * yr * tukr);
    pgrad.dufldk[3] = f * (-xr * tkz5 - zr * tkx5 - 2 * xr * zr * tukr);
    pgrad.dufldk[4] = f * (-yr * tkz5 - zr * tky5 - 2 * yr * zr * tukr);
    pgrad.dufldk[5] = f * (-zr * tkz5 - zr * zr * tukr);

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

    // get the dEd/dR terms for Thole direct polarization force

    real depx, depy, depz;

    depx = tixx * ukxp + tixy * ukyp + tixz * ukzp - tkxx * uixp - tkxy * uiyp -
        tkxz * uizp;
    depy = tixy * ukxp + tiyy * ukyp + tiyz * ukzp - tkxy * uixp - tkyy * uiyp -
        tkyz * uizp;
    depz = tixz * ukxp + tiyz * ukyp + tizz * ukzp - tkxz * uixp - tkyz * uiyp -
        tkzz * uizp;
    if_constexpr(ETYP == elec_t::ewald) {
      pgrad.frcx += depx;
      pgrad.frcy += depy;
      pgrad.frcz += depz;
    }
    else if_constexpr(ETYP == elec_t::coulomb) {
      pgrad.frcx = dscale * depx;
      pgrad.frcy = dscale * depy;
      pgrad.frcz = dscale * depz;
    }

    // get the dEp/dR terms for Thole direct polarization force

    depx = tixx * ukx + tixy * uky + tixz * ukz - tkxx * uix - tkxy * uiy -
        tkxz * uiz;
    depy = tixy * ukx + tiyy * uky + tiyz * ukz - tkxy * uix - tkyy * uiy -
        tkyz * uiz;
    depz = tixz * ukx + tiyz * uky + tizz * ukz - tkxz * uix - tkyz * uiy -
        tkzz * uiz;
    if_constexpr(ETYP == elec_t::ewald) {
      pgrad.frcx += depx;
      pgrad.frcy += depy;
      pgrad.frcz += depz;
    }
    else if_constexpr(ETYP == elec_t::coulomb) {
      pgrad.frcx += pscale * depx;
      pgrad.frcy += pscale * depy;
      pgrad.frcz += pscale * depz;
    }

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

    depx = tixx * ukxp + tixy * ukyp + tixz * ukzp + tkxx * uixp + tkxy * uiyp +
        tkxz * uizp;
    depy = tixy * ukxp + tiyy * ukyp + tiyz * ukzp + tkxy * uixp + tkyy * uiyp +
        tkyz * uizp;
    depz = tixz * ukxp + tiyz * ukyp + tizz * ukzp + tkxz * uixp + tkyz * uiyp +
        tkzz * uizp;
    if_constexpr(ETYP == elec_t::ewald) {
      pgrad.frcx += depx;
      pgrad.frcy += depy;
      pgrad.frcz += depz;
    }
    else if_constexpr(ETYP == elec_t::coulomb) {
      pgrad.frcx += uscale * depx;
      pgrad.frcy += uscale * depy;
      pgrad.frcz += uscale * depz;
    }

    pgrad.frcx *= f;
    pgrad.frcy *= f;
    pgrad.frcz *= f;
  }
}
TINKER_NAMESPACE_END

#endif
