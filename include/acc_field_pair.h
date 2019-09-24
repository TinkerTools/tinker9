#ifndef TINKER_ACC_FIELD_PAIR_H_
#define TINKER_ACC_FIELD_PAIR_H_

#include "acc_mathfunc.h"
#include "e_mpole.h"
#include "e_polar.h"
#include "md.h"

TINKER_NAMESPACE_BEGIN
struct FieldPair {
  real fid[3], fkd[3], fip[3], fkp[3];
};

#pragma acc routine seq
template <elec_t ETYP>
inline void dfield_pair_acc(                                          //
    real r2, real xr, real yr, real zr,                               //
    real dscale, real pscale, real aewald,                            //
    real ci, real dix, real diy, real diz,                            //
    real qixx, real qixy, real qixz, real qiyy, real qiyz, real qizz, //
    real pdi, real pti,                                               //
    real ck, real dkx, real dky, real dkz,                            //
    real qkxx, real qkxy, real qkxz, real qkyy, real qkyz, real qkzz, //
    real pdk, real ptk,                                               //
    FieldPair& __restrict__ pairf) {
  real r = REAL_SQRT(r2);
  real invr1 = REAL_RECIP(r);
  real rr2 = invr1 * invr1;

  real scale3 = 1;
  real scale5 = 1;
  real scale7 = 1;
  // if use_thole
  real damp = pdi * pdk;
  if (damp != 0) {
    real pgamma = REAL_MIN(pti, ptk);
    damp = -pgamma * REAL_CUBE(r / damp);
    if (damp > -50) {
      real expdamp = REAL_EXP(damp);
      scale3 = 1 - expdamp;
      scale5 = 1 - expdamp * (1 - damp);
      scale7 = 1 - expdamp * (1 - damp + (real)0.6 * REAL_SQ(damp));
    }
  }

  real bn[4], rr1, rr3, rr5, rr7;
  if_constexpr(ETYP == elec_t::ewald) {
    real ralpha = aewald * r;
    bn[0] = REAL_ERFC(ralpha) * invr1;
    real alsq2 = 2 * REAL_SQ(aewald);
    real alsq2n = (aewald > 0 ? REAL_RECIP(sqrtpi * aewald) : 0);
    real exp2a = REAL_EXP(-REAL_SQ(ralpha));
    for (int j = 1; j <= 3; ++j) {
      alsq2n *= alsq2;
      bn[j] = ((j + j - 1) * bn[j - 1] + alsq2n * exp2a) * rr2;
    }
  }
  rr1 = invr1;
  rr3 = rr1 * rr2;
  rr5 = 3 * rr1 * rr2 * rr2;
  rr7 = 15 * rr1 * rr2 * rr2 * rr2;

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

  real bcn1, bcn2, bcn3;

  // d-field

  if_constexpr(ETYP == elec_t::ewald) {
    bcn1 = bn[1] - (1 - scale3) * rr3;
    bcn2 = bn[2] - (1 - scale5) * rr5;
    bcn3 = bn[3] - (1 - scale7) * rr7;
  }
  else if_constexpr(ETYP == elec_t::coulomb) {
    bcn1 = dscale * scale3 * rr3;
    bcn2 = dscale * scale5 * rr5;
    bcn3 = dscale * scale7 * rr7;
  }

  pairf.fid[0] =
      -xr * (bcn1 * ck - bcn2 * dkr + bcn3 * qkr) - bcn1 * dkx + 2 * bcn2 * qkx;
  pairf.fid[1] =
      -yr * (bcn1 * ck - bcn2 * dkr + bcn3 * qkr) - bcn1 * dky + 2 * bcn2 * qky;
  pairf.fid[2] =
      -zr * (bcn1 * ck - bcn2 * dkr + bcn3 * qkr) - bcn1 * dkz + 2 * bcn2 * qkz;
  pairf.fkd[0] =
      xr * (bcn1 * ci + bcn2 * dir + bcn3 * qir) - bcn1 * dix - 2 * bcn2 * qix;
  pairf.fkd[1] =
      yr * (bcn1 * ci + bcn2 * dir + bcn3 * qir) - bcn1 * diy - 2 * bcn2 * qiy;
  pairf.fkd[2] =
      zr * (bcn1 * ci + bcn2 * dir + bcn3 * qir) - bcn1 * diz - 2 * bcn2 * qiz;

  // p-field

  if_constexpr(ETYP == elec_t::ewald) {
    pairf.fip[0] = pairf.fid[0];
    pairf.fip[1] = pairf.fid[1];
    pairf.fip[2] = pairf.fid[2];
    pairf.fkp[0] = pairf.fkd[0];
    pairf.fkp[1] = pairf.fkd[1];
    pairf.fkp[2] = pairf.fkd[2];
  }
  else if_constexpr(ETYP == elec_t::coulomb) {
    bcn1 = pscale * scale3 * rr3;
    bcn2 = pscale * scale5 * rr5;
    bcn3 = pscale * scale7 * rr7;
    pairf.fip[0] = -xr * (bcn1 * ck - bcn2 * dkr + bcn3 * qkr) - bcn1 * dkx +
        2 * bcn2 * qkx;
    pairf.fip[1] = -yr * (bcn1 * ck - bcn2 * dkr + bcn3 * qkr) - bcn1 * dky +
        2 * bcn2 * qky;
    pairf.fip[2] = -zr * (bcn1 * ck - bcn2 * dkr + bcn3 * qkr) - bcn1 * dkz +
        2 * bcn2 * qkz;
    pairf.fkp[0] = xr * (bcn1 * ci + bcn2 * dir + bcn3 * qir) - bcn1 * dix -
        2 * bcn2 * qix;
    pairf.fkp[1] = yr * (bcn1 * ci + bcn2 * dir + bcn3 * qir) - bcn1 * diy -
        2 * bcn2 * qiy;
    pairf.fkp[2] = zr * (bcn1 * ci + bcn2 * dir + bcn3 * qir) - bcn1 * diz -
        2 * bcn2 * qiz;
  }
}

template <int USE, elec_t ETYP>
inline void ufield_pair_acc(               //
    real r2, real xr, real yr, real zr,    //
    real dscale, real pscale, real aewald, //
    real uindi0, real uindi1, real uindi2, //
    real uinpi0, real uinpi1, real uinpi2, //
    real pdi, real pti,                    //
    real uindk0, real uindk1, real uindk2, //
    real uinpk0, real uinpk1, real uinpk2, //
    real pdk, real ptk,                    //
    FieldPair& __restrict__ pairf) {}
TINKER_NAMESPACE_END

#endif
