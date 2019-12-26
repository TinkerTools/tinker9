#include "add.h"
#include "e_mplar.h"
#include "e_mpole.h"
#include "e_polar.h"
#include "empole_self.h"
#include "epolar_trq.h"
#include "launch.h"
#include "md.h"
#include "pme.h"
#include "seq_damp.h"
#include "seq_image.h"
#include "spatial.h"


TINKER_NAMESPACE_BEGIN
template <int USE, elec_t ETYP>
__device__
void pair_mplar(
   real r2, real xr, real yr, real zr, real mscale, real dscale, real pscale,
   real uscale, //
   real ci, real dix, real diy, real diz, real qixx, real qixy, real qixz,
   real qiyy, real qiyz, real qizz, real uix, real uiy, real uiz, real uixp,
   real uiyp, real uizp, real pdi, real pti, //
   real ck, real dkx, real dky, real dkz, real qkxx, real qkxy, real qkxz,
   real qkyy, real qkyz, real qkzz, real ukx, real uky, real ukz, real ukxp,
   real ukyp, real ukzp, real pdk, real ptk, //
   real f, real aewald,                      //
   real3& restrict frci, real3& restrict frck, real3& restrict trqi,
   real3& restrict trqk, real3& restrict ufldi, real3& restrict ufldk,
   real (&restrict dufldi)[6], real (&restrict dufldk)[6], //
   real& restrict etl, real& restrict vtlxx, real& restrict vtlxy,
   real& restrict vtlxz, real& restrict vtlyy, real& restrict vtlyz,
   real& restrict vtlzz)
{
   constexpr int do_e = USE & calc::energy;
   constexpr int do_g = USE & calc::grad;
   constexpr int do_v = USE & calc::virial;


   real r = REAL_SQRT(r2);
   real invr1 = REAL_RECIP(r);
   real rr2 = invr1 * invr1;
   real rr1 = invr1;
   real rr3 = rr1 * rr2;
   real rr5 = 3 * rr3 * rr2;
   real rr7 = 5 * rr5 * rr2;
   real rr9 = 7 * rr7 * rr2;
   real rr11;
   if CONSTEXPR (do_g) {
      rr11 = 9 * rr9 * rr2;
   }
   real bn[6];
   if CONSTEXPR (ETYP == elec_t::ewald) {
      if CONSTEXPR (!do_g) {
         damp_ewald<5>(bn, r, invr1, rr2, aewald);
      } else {
         damp_ewald<6>(bn, r, invr1, rr2, aewald);
      }
   } else if CONSTEXPR (ETYP == elec_t::coulomb) {
      bn[0] = rr1;
      bn[1] = rr3;
      bn[2] = rr5;
      bn[3] = rr7;
      bn[4] = rr9;
      if CONSTEXPR (do_g) {
         bn[5] = rr11;
      }
   }


   real e;
   real3 frc;
   if CONSTEXPR (do_e)
      e = 0;
   if CONSTEXPR (do_g)
      frc = make_real3(0, 0, 0);


   // empole


   if CONSTEXPR (ETYP == elec_t::ewald) {
      mscale = 1;
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
   real dik = dix * dkx + diy * dky + diz * dkz;
   real qik = qix * qkx + qiy * qky + qiz * qkz;
   real diqk = dix * qkx + diy * qky + diz * qkz;
   real dkqi = dkx * qix + dky * qiy + dkz * qiz;
   real qiqk = 2 * (qixy * qkxy + qixz * qkxz + qiyz * qkyz) + qixx * qkxx +
      qiyy * qkyy + qizz * qkzz;


   real term1 = ci * ck;
   real term2 = ck * dir - ci * dkr + dik;
   real term3 = ci * qkr + ck * qir - dir * dkr + 2 * (dkqi - diqk + qiqk);
   real term4 = dir * qkr - dkr * qir - 4 * qik;
   real term5 = qir * qkr;


   if CONSTEXPR (do_e) {
      e += mscale * f *
         (term1 * bn[0] + term2 * bn[1] + term3 * bn[2] + term4 * bn[3] +
          term5 * bn[4]);
   }


   if CONSTEXPR (do_g) {


      // gradient


      real qixk = qixx * qkx + qixy * qky + qixz * qkz;
      real qiyk = qixy * qkx + qiyy * qky + qiyz * qkz;
      real qizk = qixz * qkx + qiyz * qky + qizz * qkz;
      real qkxi = qkxx * qix + qkxy * qiy + qkxz * qiz;
      real qkyi = qkxy * qix + qkyy * qiy + qkyz * qiz;
      real qkzi = qkxz * qix + qkyz * qiy + qkzz * qiz;

      real diqkx = dix * qkxx + diy * qkxy + diz * qkxz;
      real diqky = dix * qkxy + diy * qkyy + diz * qkyz;
      real diqkz = dix * qkxz + diy * qkyz + diz * qkzz;
      real dkqix = dkx * qixx + dky * qixy + dkz * qixz;
      real dkqiy = dkx * qixy + dky * qiyy + dkz * qiyz;
      real dkqiz = dkx * qixz + dky * qiyz + dkz * qizz;


      real de = term1 * bn[1] + term2 * bn[2] + term3 * bn[3] + term4 * bn[4] +
         term5 * bn[5];


      term1 = -ck * bn[1] + dkr * bn[2] - qkr * bn[3];
      term2 = ci * bn[1] + dir * bn[2] + qir * bn[3];
      term3 = 2 * bn[2];
      term4 = 2 * (-ck * bn[2] + dkr * bn[3] - qkr * bn[4]);
      term5 = 2 * (-ci * bn[2] - dir * bn[3] - qir * bn[4]);
      real term6 = 4 * bn[3];


      real3 frc0;
      frc0.x = de * xr + term1 * dix + term2 * dkx + term3 * (diqkx - dkqix) +
         term4 * qix + term5 * qkx + term6 * (qixk + qkxi);
      frc0.y = de * yr + term1 * diy + term2 * dky + term3 * (diqky - dkqiy) +
         term4 * qiy + term5 * qky + term6 * (qiyk + qkyi);
      frc0.z = de * zr + term1 * diz + term2 * dkz + term3 * (diqkz - dkqiz) +
         term4 * qiz + term5 * qkz + term6 * (qizk + qkzi);
      frc += (mscale * f) * frc0;


      // torque


      real dirx = diy * zr - diz * yr;
      real diry = diz * xr - dix * zr;
      real dirz = dix * yr - diy * xr;
      real dkrx = dky * zr - dkz * yr;
      real dkry = dkz * xr - dkx * zr;
      real dkrz = dkx * yr - dky * xr;
      real dikx = diy * dkz - diz * dky;
      real diky = diz * dkx - dix * dkz;
      real dikz = dix * dky - diy * dkx;

      real qirx = qiz * yr - qiy * zr;
      real qiry = qix * zr - qiz * xr;
      real qirz = qiy * xr - qix * yr;
      real qkrx = qkz * yr - qky * zr;
      real qkry = qkx * zr - qkz * xr;
      real qkrz = qky * xr - qkx * yr;
      real qikx = qky * qiz - qkz * qiy;
      real qiky = qkz * qix - qkx * qiz;
      real qikz = qkx * qiy - qky * qix;

      real qikrx = qizk * yr - qiyk * zr;
      real qikry = qixk * zr - qizk * xr;
      real qikrz = qiyk * xr - qixk * yr;
      real qkirx = qkzi * yr - qkyi * zr;
      real qkiry = qkxi * zr - qkzi * xr;
      real qkirz = qkyi * xr - qkxi * yr;

      real diqkrx = diqkz * yr - diqky * zr;
      real diqkry = diqkx * zr - diqkz * xr;
      real diqkrz = diqky * xr - diqkx * yr;
      real dkqirx = dkqiz * yr - dkqiy * zr;
      real dkqiry = dkqix * zr - dkqiz * xr;
      real dkqirz = dkqiy * xr - dkqix * yr;


      real dqikx = diy * qkz - diz * qky + dky * qiz - dkz * qiy -
         2 *
            (qixy * qkxz + qiyy * qkyz + qiyz * qkzz - qixz * qkxy -
             qiyz * qkyy - qizz * qkyz);
      real dqiky = diz * qkx - dix * qkz + dkz * qix - dkx * qiz -
         2 *
            (qixz * qkxx + qiyz * qkxy + qizz * qkxz - qixx * qkxz -
             qixy * qkyz - qixz * qkzz);
      real dqikz = dix * qky - diy * qkx + dkx * qiy - dky * qix -
         2 *
            (qixx * qkxy + qixy * qkyy + qixz * qkyz - qixy * qkxx -
             qiyy * qkxy - qiyz * qkxz);


      real3 trq0;
      trq0.x = -bn[1] * dikx + term1 * dirx + term3 * (dqikx + dkqirx) -
         term4 * qirx - term6 * (qikrx + qikx);
      trq0.y = -bn[1] * diky + term1 * diry + term3 * (dqiky + dkqiry) -
         term4 * qiry - term6 * (qikry + qiky);
      trq0.z = -bn[1] * dikz + term1 * dirz + term3 * (dqikz + dkqirz) -
         term4 * qirz - term6 * (qikrz + qikz);
      trqi += (mscale * f) * trq0;
      trq0.x = bn[1] * dikx + term2 * dkrx - term3 * (dqikx + diqkrx) -
         term5 * qkrx - term6 * (qkirx - qikx);
      trq0.y = bn[1] * diky + term2 * dkry - term3 * (dqiky + diqkry) -
         term5 * qkry - term6 * (qkiry - qiky);
      trq0.z = bn[1] * dikz + term2 * dkrz - term3 * (dqikz + diqkrz) -
         term5 * qkrz - term6 * (qkirz - qikz);
      trqk += (mscale * f) * trq0;
   }


   // epolar


   if CONSTEXPR (ETYP == elec_t::ewald) {
      dscale = 1;
      pscale = 1;
      uscale = 1;
   }


   f *= 0.5f;
   real uir = uix * xr + uiy * yr + uiz * zr;
   real ukr = ukx * xr + uky * yr + ukz * zr;


   // if use_thole
   real ex3, ex5, ex7;
   real rc31, rc32, rc33, rc51, rc52, rc53, rc71, rc72, rc73;
   if CONSTEXPR (!do_g) {
      damp_thole3(r, pdi, pti, pdk, ptk, ex3, ex5, ex7);
      ex3 = 1 - ex3;
      ex5 = 1 - ex5;
      ex7 = 1 - ex7;
   } else {
      damp_thole3g(r, rr2, xr, yr, zr, pdi, pti, pdk, ptk, ex3, ex5, ex7, rc31,
                   rc32, rc33, rc51, rc52, rc53, rc71, rc72, rc73);
      rc31 *= rr3;
      rc32 *= rr3;
      rc33 *= rr3;
      rc51 *= rr5;
      rc52 *= rr5;
      rc53 *= rr5;
      rc71 *= rr7;
      rc72 *= rr7;
      rc73 *= rr7;
   }
   // end if use_thole


   if CONSTEXPR (do_g) {
      real uirp = uixp * xr + uiyp * yr + uizp * zr;
      real ukrp = ukxp * xr + ukyp * yr + ukzp * zr;


      real sr3 = bn[1] - ex3 * rr3;
      real sr5 = bn[2] - ex5 * rr5;
      real sr7 = bn[3] - ex7 * rr7;


      // get the induced dipole field used for dipole torques


      real tuir, tukr;


      real tix3 = pscale * sr3 * ukx + dscale * sr3 * ukxp;
      real tiy3 = pscale * sr3 * uky + dscale * sr3 * ukyp;
      real tiz3 = pscale * sr3 * ukz + dscale * sr3 * ukzp;
      real tkx3 = pscale * sr3 * uix + dscale * sr3 * uixp;
      real tky3 = pscale * sr3 * uiy + dscale * sr3 * uiyp;
      real tkz3 = pscale * sr3 * uiz + dscale * sr3 * uizp;
      tuir = -pscale * sr5 * ukr - dscale * sr5 * ukrp;
      tukr = -pscale * sr5 * uir - dscale * sr5 * uirp;


      ufldi.x += f * (tix3 + xr * tuir);
      ufldi.y += f * (tiy3 + yr * tuir);
      ufldi.z += f * (tiz3 + zr * tuir);
      ufldk.x += f * (tkx3 + xr * tukr);
      ufldk.y += f * (tky3 + yr * tukr);
      ufldk.z += f * (tkz3 + zr * tukr);


      // get induced dipole field gradient used for quadrupole torques


      real tix5 = 2 * (pscale * sr5 * ukx + dscale * sr5 * ukxp);
      real tiy5 = 2 * (pscale * sr5 * uky + dscale * sr5 * ukyp);
      real tiz5 = 2 * (pscale * sr5 * ukz + dscale * sr5 * ukzp);
      real tkx5 = 2 * (pscale * sr5 * uix + dscale * sr5 * uixp);
      real tky5 = 2 * (pscale * sr5 * uiy + dscale * sr5 * uiyp);
      real tkz5 = 2 * (pscale * sr5 * uiz + dscale * sr5 * uizp);
      tuir = -pscale * sr7 * ukr - dscale * sr7 * ukrp;
      tukr = -pscale * sr7 * uir - dscale * sr7 * uirp;


      dufldi[0] += f * (xr * tix5 + xr * xr * tuir);
      dufldi[1] += f * (xr * tiy5 + yr * tix5 + 2 * xr * yr * tuir);
      dufldi[2] += f * (yr * tiy5 + yr * yr * tuir);
      dufldi[3] += f * (xr * tiz5 + zr * tix5 + 2 * xr * zr * tuir);
      dufldi[4] += f * (yr * tiz5 + zr * tiy5 + 2 * yr * zr * tuir);
      dufldi[5] += f * (zr * tiz5 + zr * zr * tuir);
      dufldk[0] += f * (-xr * tkx5 - xr * xr * tukr);
      dufldk[1] += f * (-xr * tky5 - yr * tkx5 - 2 * xr * yr * tukr);
      dufldk[2] += f * (-yr * tky5 - yr * yr * tukr);
      dufldk[3] += f * (-xr * tkz5 - zr * tkx5 - 2 * xr * zr * tukr);
      dufldk[4] += f * (-yr * tkz5 - zr * tky5 - 2 * yr * zr * tukr);
      dufldk[5] += f * (-zr * tkz5 - zr * zr * tukr);


      // get the field gradient for direct polarization force


      real term1, term2, term3, term4, term5, term6, term7;


      term1 = bn[2] - ex3 * rr5;
      term2 = bn[3] - ex5 * rr7;
      term3 = -sr3 + term1 * xr * xr - xr * rc31;
      term4 = rc31 - term1 * xr - sr5 * xr;
      term5 = term2 * xr * xr - sr5 - xr * rc51;
      term6 = (bn[4] - ex7 * rr9) * xr * xr - bn[3] - xr * rc71;
      term7 = rc51 - 2 * bn[3] * xr + (ex5 + 1.5f * ex7) * rr7 * xr;
      real tixx = ci * term3 + dix * term4 + dir * term5 + 2 * sr5 * qixx +
         (qiy * yr + qiz * zr) * ex7 * rr7 + 2 * qix * term7 + qir * term6;
      real tkxx = ck * term3 - dkx * term4 - dkr * term5 + 2 * sr5 * qkxx +
         (qky * yr + qkz * zr) * ex7 * rr7 + 2 * qkx * term7 + qkr * term6;


      term3 = -sr3 + term1 * yr * yr - yr * rc32;
      term4 = rc32 - term1 * yr - sr5 * yr;
      term5 = term2 * yr * yr - sr5 - yr * rc52;
      term6 = (bn[4] - ex7 * rr9) * yr * yr - bn[3] - yr * rc72;
      term7 = rc52 - 2 * bn[3] * yr + (ex5 + 1.5f * ex7) * rr7 * yr;
      real tiyy = ci * term3 + diy * term4 + dir * term5 + 2 * sr5 * qiyy +
         (qix * xr + qiz * zr) * ex7 * rr7 + 2 * qiy * term7 + qir * term6;
      real tkyy = ck * term3 - dky * term4 - dkr * term5 + 2 * sr5 * qkyy +
         (qkx * xr + qkz * zr) * ex7 * rr7 + 2 * qky * term7 + qkr * term6;


      term3 = -sr3 + term1 * zr * zr - zr * rc33;
      term4 = rc33 - term1 * zr - sr5 * zr;
      term5 = term2 * zr * zr - sr5 - zr * rc53;
      term6 = (bn[4] - ex7 * rr9) * zr * zr - bn[3] - zr * rc73;
      term7 = rc53 - 2 * bn[3] * zr + (ex5 + 1.5f * ex7) * rr7 * zr;
      real tizz = ci * term3 + diz * term4 + dir * term5 + 2 * sr5 * qizz +
         (qix * xr + qiy * yr) * ex7 * rr7 + 2 * qiz * term7 + qir * term6;
      real tkzz = ck * term3 - dkz * term4 - dkr * term5 + 2 * sr5 * qkzz +
         (qkx * xr + qky * yr) * ex7 * rr7 + 2 * qkz * term7 + qkr * term6;


      term3 = term1 * xr * yr - yr * rc31;
      term4 = rc31 - term1 * xr;
      term5 = term2 * xr * yr - yr * rc51;
      term6 = (bn[4] - ex7 * rr9) * xr * yr - yr * rc71;
      term7 = rc51 - term2 * xr;
      real tixy = ci * term3 - sr5 * dix * yr + diy * term4 + dir * term5 +
         2 * sr5 * qixy - 2 * sr7 * yr * qix + 2 * qiy * term7 + qir * term6;
      real tkxy = ck * term3 + sr5 * dkx * yr - dky * term4 - dkr * term5 +
         2 * sr5 * qkxy - 2 * sr7 * yr * qkx + 2 * qky * term7 + qkr * term6;


      term3 = term1 * xr * zr - zr * rc31;
      term5 = term2 * xr * zr - zr * rc51;
      term6 = (bn[4] - ex7 * rr9) * xr * zr - zr * rc71;
      real tixz = ci * term3 - sr5 * dix * zr + diz * term4 + dir * term5 +
         2 * sr5 * qixz - 2 * sr7 * zr * qix + 2 * qiz * term7 + qir * term6;
      real tkxz = ck * term3 + sr5 * dkx * zr - dkz * term4 - dkr * term5 +
         2 * sr5 * qkxz - 2 * sr7 * zr * qkx + 2 * qkz * term7 + qkr * term6;


      term3 = term1 * yr * zr - zr * rc32;
      term4 = rc32 - term1 * yr;
      term5 = term2 * yr * zr - zr * rc52;
      term6 = (bn[4] - ex7 * rr9) * yr * zr - zr * rc72;
      term7 = rc52 - term2 * yr;
      real tiyz = ci * term3 - sr5 * diy * zr + diz * term4 + dir * term5 +
         2 * sr5 * qiyz - 2 * sr7 * zr * qiy + 2 * qiz * term7 + qir * term6;
      real tkyz = ck * term3 + sr5 * dky * zr - dkz * term4 - dkr * term5 +
         2 * sr5 * qkyz - 2 * sr7 * zr * qky + 2 * qkz * term7 + qkr * term6;


      // get the dEd/dR terms for Thole direct polarization force


      real depx, depy, depz;


      depx = tixx * ukxp + tixy * ukyp + tixz * ukzp - tkxx * uixp -
         tkxy * uiyp - tkxz * uizp;
      depy = tixy * ukxp + tiyy * ukyp + tiyz * ukzp - tkxy * uixp -
         tkyy * uiyp - tkyz * uizp;
      depz = tixz * ukxp + tiyz * ukyp + tizz * ukzp - tkxz * uixp -
         tkyz * uiyp - tkzz * uizp;
      if CONSTEXPR (ETYP == elec_t::ewald) {
         frc.x += f * -depx;
         frc.y += f * -depy;
         frc.z += f * -depz;
      } else if CONSTEXPR (ETYP == elec_t::coulomb) {
         frc.x += f * -depx * dscale;
         frc.y += f * -depy * dscale;
         frc.z += f * -depz * dscale;
      }


      // get the dEp/dR terms for Thole direct polarization force


      depx = tixx * ukx + tixy * uky + tixz * ukz - tkxx * uix - tkxy * uiy -
         tkxz * uiz;
      depy = tixy * ukx + tiyy * uky + tiyz * ukz - tkxy * uix - tkyy * uiy -
         tkyz * uiz;
      depz = tixz * ukx + tiyz * uky + tizz * ukz - tkxz * uix - tkyz * uiy -
         tkzz * uiz;
      if CONSTEXPR (ETYP == elec_t::ewald) {
         frc.x -= f * depx;
         frc.y -= f * depy;
         frc.z -= f * depz;
      } else if CONSTEXPR (ETYP == elec_t::coulomb) {
         frc.x -= f * pscale * depx;
         frc.y -= f * pscale * depy;
         frc.z -= f * pscale * depz;
      }


      // get the dtau/dr terms used for mutual polarization force


      term1 = bn[2] - ex3 * rr5;
      term2 = bn[3] - ex5 * rr7;
      term3 = sr5 + term1;


      term5 = -xr * term3 + rc31;
      term6 = -sr5 + xr * xr * term2 - xr * rc51;
      tixx = uix * term5 + uir * term6;
      tkxx = ukx * term5 + ukr * term6;


      term5 = -yr * term3 + rc32;
      term6 = -sr5 + yr * yr * term2 - yr * rc52;
      tiyy = uiy * term5 + uir * term6;
      tkyy = uky * term5 + ukr * term6;


      term5 = -zr * term3 + rc33;
      term6 = -sr5 + zr * zr * term2 - zr * rc53;
      tizz = uiz * term5 + uir * term6;
      tkzz = ukz * term5 + ukr * term6;


      term4 = -sr5 * yr;
      term5 = -xr * term1 + rc31;
      term6 = xr * yr * term2 - yr * rc51;
      tixy = uix * term4 + uiy * term5 + uir * term6;
      tkxy = ukx * term4 + uky * term5 + ukr * term6;


      term4 = -sr5 * zr;
      term6 = xr * zr * term2 - zr * rc51;
      tixz = uix * term4 + uiz * term5 + uir * term6;
      tkxz = ukx * term4 + ukz * term5 + ukr * term6;


      term5 = -yr * term1 + rc32;
      term6 = yr * zr * term2 - zr * rc52;
      tiyz = uiy * term4 + uiz * term5 + uir * term6;
      tkyz = uky * term4 + ukz * term5 + ukr * term6;


      depx = tixx * ukxp + tixy * ukyp + tixz * ukzp + tkxx * uixp +
         tkxy * uiyp + tkxz * uizp;
      depy = tixy * ukxp + tiyy * ukyp + tiyz * ukzp + tkxy * uixp +
         tkyy * uiyp + tkyz * uizp;
      depz = tixz * ukxp + tiyz * ukyp + tizz * ukzp + tkxz * uixp +
         tkyz * uiyp + tkzz * uizp;
      if CONSTEXPR (ETYP == elec_t::ewald) {
         frc.x -= f * depx;
         frc.y -= f * depy;
         frc.z -= f * depz;
      } else if CONSTEXPR (ETYP == elec_t::coulomb) {
         frc.x -= f * uscale * depx;
         frc.y -= f * uscale * depy;
         frc.z -= f * uscale * depz;
      }
   }


   // save results


   if CONSTEXPR (do_e) {
      etl += e;
   }
   if CONSTEXPR (do_g) {
      frci += frc;
      frck -= frc;
   }
   if CONSTEXPR (do_v) {
      vtlxx += -xr * frc.x;
      vtlxy += -0.5f * (yr * frc.x + xr * frc.y);
      vtlxz += -0.5f * (zr * frc.x + xr * frc.z);
      vtlyy += -yr * frc.y;
      vtlyz += -0.5f * (zr * frc.y + yr * frc.z);
      vtlzz += -zr * frc.z;
   }
}


#define EMPLAR_ARGS                                                            \
   size_t bufsize, energy_buffer restrict ebuf,                                \
      virial_buffer restrict vir_ebuf, real *restrict gx, real *restrict gy,   \
      real *restrict gz, real *restrict trqx, real *restrict trqy,             \
      real *restrict trqz, TINKER_IMAGE_PARAMS, real off, real f,              \
      const real(*restrict rpole)[10], const real *restrict pdamp,             \
      const real *restrict thole, const real(*restrict uind)[3],               \
      const real(*restrict uinp)[3], real(*restrict ufld)[3],                  \
      real(*restrict dufld)[6]


template <int USE, elec_t ETYP>
__launch_bounds__(BLOCK_DIM) __global__
void emplar_cu1(EMPLAR_ARGS, int n, const Spatial::SortedAtom* restrict sorted,
                int niak, const int* restrict iak, const int* restrict lst,
                real aewald)
{
   constexpr int do_e = USE & calc::energy;
   constexpr int do_g = USE & calc::grad;
   constexpr int do_v = USE & calc::virial;


   const int iwarp = (threadIdx.x + blockIdx.x * blockDim.x) / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);
   const int offset = (threadIdx.x + blockIdx.x * blockDim.x) & (bufsize - 1);


   struct Data
   {
      real3 frc, trq; // force and torque
      real3 ufld;
      real dufld[6];
      real3 pos;
      real c;  // charge
      real3 d; // dipole
      real qxx, qxy, qxz, qyy, qyz, qzz;
      real3 ud, up;
      real damp, thole;
   };
   __shared__ Data data[BLOCK_DIM];


   const real off2 = off * off;
   for (int iw = iwarp; iw < niak; iw += nwarp) {
      real etl;
      if CONSTEXPR (do_e) {
         etl = 0;
      }
      real vtlxx, vtlxy, vtlxz, vtlyy, vtlyz, vtlzz;
      if CONSTEXPR (do_v) {
         vtlxx = 0;
         vtlxy = 0;
         vtlxz = 0;
         vtlyy = 0;
         vtlyz = 0;
         vtlzz = 0;
      }


      Data idat;
      if CONSTEXPR (do_g) {
         idat.frc = make_real3(0, 0, 0);
         idat.trq = make_real3(0, 0, 0);
         idat.ufld = make_real3(0, 0, 0);
         #pragma unroll
         for (int i = 0; i < 6; ++i) {
            idat.dufld[i] = 0;
         }
      }
      int atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
      idat.pos = make_real3(sorted[atomi].x, sorted[atomi].y, sorted[atomi].z);
      int i = sorted[atomi].unsorted;
      idat.c = rpole[i][mpl_pme_0];
      idat.d = make_real3(rpole[i][mpl_pme_x], rpole[i][mpl_pme_y],
                          rpole[i][mpl_pme_z]);
      idat.qxx = rpole[i][mpl_pme_xx];
      idat.qxy = rpole[i][mpl_pme_xy];
      idat.qxz = rpole[i][mpl_pme_xz];
      idat.qyy = rpole[i][mpl_pme_yy];
      idat.qyz = rpole[i][mpl_pme_yz];
      idat.qzz = rpole[i][mpl_pme_zz];
      idat.ud = make_real3(uind[i][0], uind[i][1], uind[i][2]);
      idat.up = make_real3(uinp[i][0], uinp[i][1], uinp[i][2]);
      idat.damp = pdamp[i];
      idat.thole = thole[i];


      if CONSTEXPR (do_g) {
         data[threadIdx.x].frc = make_real3(0, 0, 0);
         data[threadIdx.x].trq = make_real3(0, 0, 0);
         data[threadIdx.x].ufld = make_real3(0, 0, 0);
         #pragma unroll
         for (int i = 0; i < 6; ++i) {
            data[threadIdx.x].dufld[i] = 0;
         }
      }
      int shatomk = lst[iw * WARP_SIZE + ilane];
      data[threadIdx.x].pos =
         make_real3(sorted[shatomk].x, sorted[shatomk].y, sorted[shatomk].z);
      int shk = sorted[shatomk].unsorted;
      data[threadIdx.x].c = rpole[shk][mpl_pme_0];
      data[threadIdx.x].d = make_real3(
         rpole[shk][mpl_pme_x], rpole[shk][mpl_pme_y], rpole[shk][mpl_pme_z]);
      data[threadIdx.x].qxx = rpole[shk][mpl_pme_xx];
      data[threadIdx.x].qxy = rpole[shk][mpl_pme_xy];
      data[threadIdx.x].qxz = rpole[shk][mpl_pme_xz];
      data[threadIdx.x].qyy = rpole[shk][mpl_pme_yy];
      data[threadIdx.x].qyz = rpole[shk][mpl_pme_yz];
      data[threadIdx.x].qzz = rpole[shk][mpl_pme_zz];
      data[threadIdx.x].ud =
         make_real3(uind[shk][0], uind[shk][1], uind[shk][2]);
      data[threadIdx.x].up =
         make_real3(uinp[shk][0], uinp[shk][1], uinp[shk][2]);
      data[threadIdx.x].damp = pdamp[shk];
      data[threadIdx.x].thole = thole[shk];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
         real3 dr = data[klane].pos - idat.pos;


         real r2 = image2(dr.x, dr.y, dr.z);
         if (atomi < atomk && r2 <= off2) {
            if CONSTEXPR (ETYP == elec_t::ewald) {
               pair_mplar<USE, elec_t::ewald>(
                  r2, dr.x, dr.y, dr.z, 1, 1, 1, 1, //
                  idat.c, idat.d.x, idat.d.y, idat.d.z, idat.qxx, idat.qxy,
                  idat.qxz, idat.qyy, idat.qyz, idat.qzz, idat.ud.x, idat.ud.y,
                  idat.ud.z, idat.up.x, idat.up.y, idat.up.z, idat.damp,
                  idat.thole, //
                  data[klane].c, data[klane].d.x, data[klane].d.y,
                  data[klane].d.z, data[klane].qxx, data[klane].qxy,
                  data[klane].qxz, data[klane].qyy, data[klane].qyz,
                  data[klane].qzz, data[klane].ud.x, data[klane].ud.y,
                  data[klane].ud.z, data[klane].up.x, data[klane].up.y,
                  data[klane].up.z, data[klane].damp, data[klane].thole, //
                  f, aewald,                                             //
                  idat.frc, data[klane].frc, idat.trq, data[klane].trq,
                  idat.ufld, data[klane].ufld, idat.dufld, data[klane].dufld,
                  etl, vtlxx, vtlxy, vtlxz, vtlyy, vtlyz, vtlzz);
            }
            if CONSTEXPR (ETYP == elec_t::coulomb) {
               pair_mplar<USE, elec_t::coulomb>(
                  r2, dr.x, dr.y, dr.z, 1, 1, 1, 1, //
                  idat.c, idat.d.x, idat.d.y, idat.d.z, idat.qxx, idat.qxy,
                  idat.qxz, idat.qyy, idat.qyz, idat.qzz, idat.ud.x, idat.ud.y,
                  idat.ud.z, idat.up.x, idat.up.y, idat.up.z, idat.damp,
                  idat.thole, //
                  data[klane].c, data[klane].d.x, data[klane].d.y,
                  data[klane].d.z, data[klane].qxx, data[klane].qxy,
                  data[klane].qxz, data[klane].qyy, data[klane].qyz,
                  data[klane].qzz, data[klane].ud.x, data[klane].ud.y,
                  data[klane].ud.z, data[klane].up.x, data[klane].up.y,
                  data[klane].up.z, data[klane].damp, data[klane].thole, //
                  f, 0,                                                  //
                  idat.frc, data[klane].frc, idat.trq, data[klane].trq,
                  idat.ufld, data[klane].ufld, idat.dufld, data[klane].dufld,
                  etl, vtlxx, vtlxy, vtlxz, vtlyy, vtlyz, vtlzz);
            }
         } // end if (include)
      }


      if CONSTEXPR (do_e)
         atomic_add(etl, ebuf, offset);
      if CONSTEXPR (do_g) {
         atomic_add(idat.frc.x, &gx[i]);
         atomic_add(idat.frc.y, &gy[i]);
         atomic_add(idat.frc.z, &gz[i]);
         atomic_add(data[threadIdx.x].frc.x, &gx[shk]);
         atomic_add(data[threadIdx.x].frc.y, &gy[shk]);
         atomic_add(data[threadIdx.x].frc.z, &gz[shk]);
         atomic_add(idat.trq.x, &trqx[i]);
         atomic_add(idat.trq.y, &trqy[i]);
         atomic_add(idat.trq.z, &trqz[i]);
         atomic_add(data[threadIdx.x].trq.x, &trqx[shk]);
         atomic_add(data[threadIdx.x].trq.y, &trqy[shk]);
         atomic_add(data[threadIdx.x].trq.z, &trqz[shk]);
         atomic_add(idat.ufld.x, &ufld[i][0]);
         atomic_add(idat.ufld.y, &ufld[i][1]);
         atomic_add(idat.ufld.z, &ufld[i][2]);
         atomic_add(data[threadIdx.x].ufld.x, &ufld[shk][0]);
         atomic_add(data[threadIdx.x].ufld.y, &ufld[shk][1]);
         atomic_add(data[threadIdx.x].ufld.z, &ufld[shk][2]);
         atomic_add(idat.dufld[0], &dufld[i][0]);
         atomic_add(idat.dufld[1], &dufld[i][1]);
         atomic_add(idat.dufld[2], &dufld[i][2]);
         atomic_add(idat.dufld[3], &dufld[i][3]);
         atomic_add(idat.dufld[4], &dufld[i][4]);
         atomic_add(idat.dufld[5], &dufld[i][5]);
         atomic_add(data[threadIdx.x].dufld[0], &dufld[shk][0]);
         atomic_add(data[threadIdx.x].dufld[1], &dufld[shk][1]);
         atomic_add(data[threadIdx.x].dufld[2], &dufld[shk][2]);
         atomic_add(data[threadIdx.x].dufld[3], &dufld[shk][3]);
         atomic_add(data[threadIdx.x].dufld[4], &dufld[shk][4]);
         atomic_add(data[threadIdx.x].dufld[5], &dufld[shk][5]);
      }
      if CONSTEXPR (do_v)
         atomic_add(vtlxx, vtlxy, vtlxz, vtlyy, vtlyz, vtlzz, vir_ebuf, offset);
   } // end for (iw)
}


template <int USE>
__global__
void emplar_cu2(EMPLAR_ARGS, const real* restrict x, const real* restrict y,
                const real* restrict z, int nmdpuexclude,
                const int (*restrict mdpuexclude)[2],
                const real (*restrict mdpuexclude_scale)[4])
{
   constexpr int do_e = USE & calc::energy;
   constexpr int do_g = USE & calc::grad;
   constexpr int do_v = USE & calc::virial;


   const real off2 = off * off;
   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < nmdpuexclude;
        ii += blockDim.x * gridDim.x) {
      int offset = ii & (bufsize - 1);


      int i = mdpuexclude[ii][0];
      int k = mdpuexclude[ii][1];
      real mscale = mdpuexclude_scale[ii][0];
      real dscale = mdpuexclude_scale[ii][1];
      real pscale = mdpuexclude_scale[ii][2];
      real uscale = mdpuexclude_scale[ii][3];


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
      real uixp = uinp[i][0];
      real uiyp = uinp[i][1];
      real uizp = uinp[i][2];
      real pdi = pdamp[i];
      real pti = thole[i];


      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off2) {
         real e;
         if CONSTEXPR (do_e) {
            e = 0;
         }
         real vxx, vxy, vxz, vyy, vyz, vzz;
         if CONSTEXPR (do_v) {
            vxx = 0;
            vxy = 0;
            vxz = 0;
            vyy = 0;
            vyz = 0;
            vzz = 0;
         }
         real3 frci, frck, trqi, trqk, ufldi, ufldk;
         real dufldi[6], dufldk[6];
         if CONSTEXPR (do_g) {
            frci = make_real3(0, 0, 0);
            frck = make_real3(0, 0, 0);
            trqi = make_real3(0, 0, 0);
            trqk = make_real3(0, 0, 0);
            ufldi = make_real3(0, 0, 0);
            ufldk = make_real3(0, 0, 0);
            #pragma unroll
            for (int ia = 0; ia < 6; ++ia) {
               dufldi[ia] = 0;
            }
            #pragma unroll
            for (int ia = 0; ia < 6; ++ia) {
               dufldk[ia] = 0;
            }
         }


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
         real ukxp = uinp[k][0];
         real ukyp = uinp[k][1];
         real ukzp = uinp[k][2];
         real pdk = pdamp[k];
         real ptk = thole[k];


         pair_mplar<USE, elec_t::coulomb>(
            r2, xr, yr, zr, mscale, dscale, pscale, uscale, //
            ci, dix, diy, diz, qixx, qixy, qixz, qiyy, qiyz, qizz, uix, uiy,
            uiz, uixp, uiyp, uizp, pdi, pti, //
            ck, dkx, dky, dkz, qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, ukx, uky,
            ukz, ukxp, ukyp, ukzp, pdk, ptk,                      //
            f, 0,                                                 //
            frci, frck, trqi, trqk, ufldi, ufldk, dufldi, dufldk, //
            e, vxx, vxy, vxz, vyy, vyz, vzz);


         if CONSTEXPR (do_e)
            atomic_add(e, ebuf, offset);
         if CONSTEXPR (do_g) {
            atomic_add(frci.x, &gx[i]);
            atomic_add(frci.y, &gy[i]);
            atomic_add(frci.z, &gz[i]);
            atomic_add(frck.x, &gx[k]);
            atomic_add(frck.y, &gy[k]);
            atomic_add(frck.z, &gz[k]);
            atomic_add(trqi.x, &trqx[i]);
            atomic_add(trqi.y, &trqy[i]);
            atomic_add(trqi.z, &trqz[i]);
            atomic_add(trqk.x, &trqx[k]);
            atomic_add(trqk.y, &trqy[k]);
            atomic_add(trqk.z, &trqz[k]);
            atomic_add(ufldi.x, &ufld[i][0]);
            atomic_add(ufldi.y, &ufld[i][1]);
            atomic_add(ufldi.z, &ufld[i][2]);
            atomic_add(ufldk.x, &ufld[k][0]);
            atomic_add(ufldk.y, &ufld[k][1]);
            atomic_add(ufldk.z, &ufld[k][2]);
            atomic_add(dufldi[0], &dufld[i][0]);
            atomic_add(dufldi[1], &dufld[i][1]);
            atomic_add(dufldi[2], &dufld[i][2]);
            atomic_add(dufldi[3], &dufld[i][3]);
            atomic_add(dufldi[4], &dufld[i][4]);
            atomic_add(dufldi[5], &dufld[i][5]);
            atomic_add(dufldk[0], &dufld[k][0]);
            atomic_add(dufldk[1], &dufld[k][1]);
            atomic_add(dufldk[2], &dufld[k][2]);
            atomic_add(dufldk[3], &dufld[k][3]);
            atomic_add(dufldk[4], &dufld[k][4]);
            atomic_add(dufldk[5], &dufld[k][5]);
         }
         if CONSTEXPR (do_v)
            atomic_add(vxx, vxy, vxz, vyy, vyz, vzz, vir_ebuf, offset);
      } // end if (include)
   }
}


template <int USE, elec_t ETYP>
void emplar_tmpl_cu(const real (*uind)[3], const real (*uinp)[3])
{
   const auto& st = *mspatial_unit;
   const real off = st.cutoff;
   auto bufsize = buffer_size();


   const real f = electric / dielec;
   real aewald = 0;
   if CONSTEXPR (ETYP == elec_t::ewald) {
      assert(epme_unit == ppme_unit);
      PMEUnit pu = epme_unit;
      aewald = pu->aewald;


      if CONSTEXPR (USE & calc::energy) {
         auto ker0 = empole_self_cu<calc::energy>;
         launch_k1s(nonblk, n, ker0, //
                    bufsize, nullptr, em, rpole, n, f, aewald);
      }
   }
   if (st.niak > 0) {
      auto ker1 = emplar_cu1<USE, ETYP>;
      launch_k1s(nonblk, WARP_SIZE * st.niak, ker1, //
                 bufsize, em, vir_em, gx, gy, gz, trqx, trqy, trqz,
                 TINKER_IMAGE_ARGS, off, f, rpole, pdamp, thole, uind, uinp,
                 ufld, dufld, //
                 n, st.sorted, st.niak, st.iak, st.lst, aewald);
   }
   if (nmdpuexclude > 0) {
      auto ker2 = emplar_cu2<USE>;
      launch_k1s(nonblk, nmdpuexclude, ker2, //
                 bufsize, em, vir_em, gx, gy, gz, trqx, trqy, trqz,
                 TINKER_IMAGE_ARGS, off, f, rpole, pdamp, thole, uind, uinp,
                 ufld, dufld, //
                 x, y, z, nmdpuexclude, mdpuexclude, mdpuexclude_scale);
   }
}


template <int USE>
void empole_recip_tmpl();
extern template void empole_recip_tmpl<calc::v0>();
extern template void empole_recip_tmpl<calc::v1>();
extern template void empole_recip_tmpl<calc::v3>();
extern template void empole_recip_tmpl<calc::v4>();
extern template void empole_recip_tmpl<calc::v5>();
extern template void empole_recip_tmpl<calc::v6>();


template <int USE>
void epolar_recip_self_tmpl(const real (*)[3], const real (*)[3]);
extern template void epolar_recip_self_tmpl<calc::v0>(const real (*)[3],
                                                      const real (*)[3]);
extern template void epolar_recip_self_tmpl<calc::v1>(const real (*)[3],
                                                      const real (*)[3]);
extern template void epolar_recip_self_tmpl<calc::v3>(const real (*)[3],
                                                      const real (*)[3]);
extern template void epolar_recip_self_tmpl<calc::v4>(const real (*)[3],
                                                      const real (*)[3]);
extern template void epolar_recip_self_tmpl<calc::v5>(const real (*)[3],
                                                      const real (*)[3]);
extern template void epolar_recip_self_tmpl<calc::v6>(const real (*)[3],
                                                      const real (*)[3]);


template <int USE>
void emplar_ewald_tmpl()
{
   constexpr int do_a = USE & calc::analyz;
   constexpr int do_e = USE & calc::energy;
   constexpr int do_g = USE & calc::grad;
   constexpr int do_v = USE & calc::virial;
   static_assert(do_v ? do_g : true, "");
   static_assert(!do_a, "");


   // induce
   induce(uind, uinp);


   if CONSTEXPR (do_g) {
      device_array::zero_async(n, ufld, dufld);
   }


   // empole real self
   // epolar real gradient
   emplar_tmpl_cu<USE, elec_t::ewald>(uind, uinp);
   // epolar torque
   if CONSTEXPR (USE & calc::grad) {
      launch_k1s(nonblk, n, epolar_trq_cu, //
                 trqx, trqy, trqz, n, rpole, ufld, dufld);
   }
   if CONSTEXPR (do_e) {
      epolar0_dotprod(uind, udirp);
   }


   // empole recip
   empole_recip_tmpl<USE>();
   // epolar recip self
   epolar_recip_self_tmpl<USE>(uind, uinp);
}


template <int USE>
void emplar_coulomb_tmpl()
{
   constexpr int do_a = USE & calc::analyz;
   constexpr int do_e = USE & calc::energy;
   constexpr int do_g = USE & calc::grad;
   constexpr int do_v = USE & calc::virial;
   static_assert(do_v ? do_g : true, "");
   static_assert(!do_a, "");


   // induce
   induce(uind, uinp);


   if CONSTEXPR (do_g) {
      device_array::zero(n, ufld, dufld);
   }


   // empole and epolar
   emplar_tmpl_cu<USE, elec_t::coulomb>(uind, uinp);
   // epolar torque
   if CONSTEXPR (USE & calc::grad) {
      launch_k1(n, epolar_trq_cu, //
                trqx, trqy, trqz, n, rpole, ufld, dufld);
   }
   if CONSTEXPR (do_e) {
      epolar0_dotprod(uind, uinp);
   }
}


void emplar_cu(int vers)
{
   assert(empole_electyp == epolar_electyp);
   if (empole_electyp == elec_t::coulomb) {
      if (vers == calc::v0)
         emplar_coulomb_tmpl<calc::v0>();
      else if (vers == calc::v1)
         emplar_coulomb_tmpl<calc::v1>();
      // else if (vers == calc::v3)
      //    emplar_coulomb_tmpl<calc::v3>();
      else if (vers == calc::v4)
         emplar_coulomb_tmpl<calc::v4>();
      else if (vers == calc::v5)
         emplar_coulomb_tmpl<calc::v5>();
      else if (vers == calc::v6)
         emplar_coulomb_tmpl<calc::v6>();
   } else if (empole_electyp == elec_t::ewald) {
      if (vers == calc::v0)
         emplar_ewald_tmpl<calc::v0>();
      else if (vers == calc::v1)
         emplar_ewald_tmpl<calc::v1>();
      // else if (vers == calc::v3)
      //    emplar_ewald_tmpl<calc::v3>();
      else if (vers == calc::v4)
         emplar_ewald_tmpl<calc::v4>();
      else if (vers == calc::v5)
         emplar_ewald_tmpl<calc::v5>();
      else if (vers == calc::v6)
         emplar_ewald_tmpl<calc::v6>();
   }
}
TINKER_NAMESPACE_END
