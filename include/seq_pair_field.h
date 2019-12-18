#pragma once
#include "elec.h"
#include "macro_void_cuda_def.h"
#include "md.h"
#include "seq_damp.h"


TINKER_NAMESPACE_BEGIN
struct PairField
{
   real fid[3], fkd[3], fip[3], fkp[3];
};


__device__
inline void zero(PairField& pairf)
{
   pairf.fid[0] = 0;
   pairf.fid[1] = 0;
   pairf.fid[2] = 0;
   pairf.fkd[0] = 0;
   pairf.fkd[1] = 0;
   pairf.fkd[2] = 0;
   pairf.fip[0] = 0;
   pairf.fip[1] = 0;
   pairf.fip[2] = 0;
   pairf.fkp[0] = 0;
   pairf.fkp[1] = 0;
   pairf.fkp[2] = 0;
}


#pragma acc routine seq
template <elec_t ETYP>
__device__
void pair_dfield(                                                //
   real r2, real xr, real yr, real zr, real dscale, real pscale, //
   real ci, real dix, real diy, real diz, real qixx, real qixy, real qixz,
   real qiyy, real qiyz, real qizz, real pdi, real pti, //
   real ck, real dkx, real dky, real dkz, real qkxx, real qkxy, real qkxz,
   real qkyy, real qkyz, real qkzz, real pdk, real ptk, //
   real aewald, PairField& restrict pairf)
{
   real r = REAL_SQRT(r2);
   real invr1 = REAL_RECIP(r);
   real rr2 = invr1 * invr1;

   real scale3, scale5, scale7;
   damp_thole3(r, pdi, pti, pdk, ptk, scale3, scale5, scale7);

   MAYBE_UNUSED real bn[4];
   if CONSTEXPR (ETYP == elec_t::ewald)
      damp_ewald<4>(bn, r, invr1, rr2, aewald);
   real rr1 = invr1;
   real rr3 = rr1 * rr2;
   real rr5 = 3 * rr1 * rr2 * rr2;
   real rr7 = 15 * rr1 * rr2 * rr2 * rr2;

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

   if CONSTEXPR (ETYP == elec_t::ewald) {
      bcn1 = bn[1] - (1 - scale3) * rr3;
      bcn2 = bn[2] - (1 - scale5) * rr5;
      bcn3 = bn[3] - (1 - scale7) * rr7;
   } else if CONSTEXPR (ETYP == elec_t::coulomb) {
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

   if CONSTEXPR (ETYP == elec_t::ewald) {
      pairf.fip[0] = pairf.fid[0];
      pairf.fip[1] = pairf.fid[1];
      pairf.fip[2] = pairf.fid[2];
      pairf.fkp[0] = pairf.fkd[0];
      pairf.fkp[1] = pairf.fkd[1];
      pairf.fkp[2] = pairf.fkd[2];
   } else if CONSTEXPR (ETYP == elec_t::coulomb) {
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


template <elec_t ETYP>
__device__
void pair_ufield(                                   //
   real r2, real xr, real yr, real zr, real uscale, //
   real uindi0, real uindi1, real uindi2, real uinpi0, real uinpi1, real uinpi2,
   real pdi, real pti, //
   real uindk0, real uindk1, real uindk2, real uinpk0, real uinpk1, real uinpk2,
   real pdk, real ptk, //
   real aewald, real3& restrict fid, real3& restrict fip, real3& restrict fkd,
   real3& restrict fkp)
{
   real r = REAL_SQRT(r2);
   real invr1 = REAL_RECIP(r);
   real rr2 = invr1 * invr1;

   real scale3, scale5;
   damp_thole2(r, pdi, pti, pdk, ptk, scale3, scale5);

   real bn[3];
   if CONSTEXPR (ETYP == elec_t::ewald)
      damp_ewald<3>(bn, r, invr1, rr2, aewald);
   real rr1 = invr1;
   real rr3 = rr1 * rr2;
   real rr5 = 3 * rr1 * rr2 * rr2;

   if CONSTEXPR (ETYP == elec_t::ewald) {
      bn[1] -= (1 - scale3) * rr3;
      bn[2] -= (1 - scale5) * rr5;
   } else if CONSTEXPR (ETYP == elec_t::coulomb) {
      bn[1] = uscale * scale3 * rr3;
      bn[2] = uscale * scale5 * rr5;
   }

   real coef;
   real3 dr = make_real3(xr, yr, zr);
   real3 uid = make_real3(uindi0, uindi1, uindi2);
   real3 uip = make_real3(uinpi0, uinpi1, uinpi2);
   real3 ukd = make_real3(uindk0, uindk1, uindk2);
   real3 ukp = make_real3(uinpk0, uinpk1, uinpk2);

   coef = bn[2] * dot3(dr, ukd);
   fid += coef * dr - bn[1] * ukd;

   coef = bn[2] * dot3(dr, ukp);
   fip += coef * dr - bn[1] * ukp;

   coef = bn[2] * dot3(dr, uid);
   fkd += coef * dr - bn[1] * uid;

   coef = bn[2] * dot3(dr, uip);
   fkp += coef * dr - bn[1] * uip;
}
TINKER_NAMESPACE_END
