#include "add.h"
#include "emplar.h"
#include "empole.h"
#include "empole_self.h"
#include "epolar.h"
#include "epolar_trq.h"
#include "image.h"
#include "launch.h"
#include "md.h"
#include "named_struct.h"
#include "pme.h"
#include "seq_damp.h"
#include "spatial.h"
#include "switch.h"


TINKER_NAMESPACE_BEGIN
template <class Ver, class ETYP>
__device__
void pair_mplar_v1(                                                       //
   real r2, real3 dr, real mscale, real dscale, real pscale, real uscale, //
   real ci, real3 di, real qixx, real qixy, real qixz, real qiyy, real qiyz,
   real qizz, real3 uid, real3 uip, real pdi, real pti, //
   real ck, real3 dk, real qkxx, real qkxy, real qkxz, real qkyy, real qkyz,
   real qkzz, real3 ukd, real3 ukp, real pdk, real ptk, //
   real f, real aewald,                                 //
   real3& restrict frci, real3& restrict frck, real3& restrict trqi,
   real3& restrict trqk, real& restrict etl, real& restrict vtlxx,
   real& restrict vtlxy, real& restrict vtlxz, real& restrict vtlyy,
   real& restrict vtlyz, real& restrict vtlzz)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;

   real bn[6];
   real sr3, sr5, sr7, sr9;
   {
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

      if CONSTEXPR (eq<ETYP, EWALD>()) {
         if CONSTEXPR (!do_g) {
            damp_ewald<5>(bn, r, invr1, rr2, aewald);
         } else {
            damp_ewald<6>(bn, r, invr1, rr2, aewald);
         }
      } else if CONSTEXPR (eq<ETYP, NON_EWALD>()) {
         bn[0] = rr1;
         bn[1] = rr3;
         bn[2] = rr5;
         bn[3] = rr7;
         bn[4] = rr9;
         if CONSTEXPR (do_g) {
            bn[5] = rr11;
         }
      }

      // if use_thole
      real ex3, ex5, ex7, ex9;
      damp_thole4(r, pdi, pti, pdk, ptk, ex3, ex5, ex7, ex9);
      sr3 = bn[1] - ex3 * rr3;
      sr5 = bn[2] - ex5 * rr5;
      sr7 = bn[3] - ex7 * rr7;
      sr9 = bn[4] - ex9 * rr9;
      // end if use_thole
   }

   real3 frc, trq1, trq2;
   if CONSTEXPR (do_g) {
      frc = make_real3(0, 0, 0);
      trq1 = make_real3(0, 0, 0);
      trq2 = make_real3(0, 0, 0);
   }

   real dir = dot3(di, dr);
   real3 qi_dr = matvec(qixx, qixy, qixz, qiyy, qiyz, qizz, dr);
   real qir = dot3(dr, qi_dr);
   real dkr = dot3(dk, dr);
   real3 qk_dr = matvec(qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, dr);
   real qkr = dot3(dr, qk_dr);

   // empole
   {
      real dik = dot3(di, dk);
      real qik = dot3(qi_dr, qk_dr);
      real diqk = dot3(di, qk_dr);
      real dkqi = dot3(dk, qi_dr);
      real qiqk = dot3(qixx, qiyy, qizz, qkxx, qkyy, qkzz) +
         2 * dot3(qixy, qixz, qiyz, qkxy, qkxz, qkyz);

      real term1 = ci * ck;
      real term2 = ck * dir - ci * dkr + dik;
      real term3 = ci * qkr + ck * qir - dir * dkr + 2 * (dkqi - diqk + qiqk);
      real term4 = dir * qkr - dkr * qir - 4 * qik;
      real term5 = qir * qkr;

      // energy
      if CONSTEXPR (do_e) {
         real e = term1 * bn[0] + term2 * bn[1] + term3 * bn[2] +
            term4 * bn[3] + term5 * bn[4];
         etl += f * mscale * e;
      }

      // gradient
      real qix = qi_dr.x;
      real qiy = qi_dr.y;
      real qiz = qi_dr.z;
      real qkx = qk_dr.x;
      real qky = qk_dr.y;
      real qkz = qk_dr.z;

      real qixk = qixx * qkx + qixy * qky + qixz * qkz; // |Qi Qk r>
      real qiyk = qixy * qkx + qiyy * qky + qiyz * qkz;
      real qizk = qixz * qkx + qiyz * qky + qizz * qkz;
      real qkxi = qkxx * qix + qkxy * qiy + qkxz * qiz; // |Qk Qi r>
      real qkyi = qkxy * qix + qkyy * qiy + qkyz * qiz;
      real qkzi = qkxz * qix + qkyz * qiy + qkzz * qiz;

      real diqkx = di.x * qkxx + di.y * qkxy + di.z * qkxz; // |Qk Di>
      real diqky = di.x * qkxy + di.y * qkyy + di.z * qkyz;
      real diqkz = di.x * qkxz + di.y * qkyz + di.z * qkzz;
      real dkqix = dk.x * qixx + dk.y * qixy + dk.z * qixz; // |Qi Dk>
      real dkqiy = dk.x * qixy + dk.y * qiyy + dk.z * qiyz;
      real dkqiz = dk.x * qixz + dk.y * qiyz + dk.z * qizz;

      real de = term1 * bn[1] + term2 * bn[2] + term3 * bn[3] + term4 * bn[4] +
         term5 * bn[5];

      term1 = -ck * bn[1] + dkr * bn[2] - qkr * bn[3];
      term2 = ci * bn[1] + dir * bn[2] + qir * bn[3];
      term3 = 2 * bn[2];
      term4 = 2 * (-ck * bn[2] + dkr * bn[3] - qkr * bn[4]);
      term5 = 2 * (-ci * bn[2] - dir * bn[3] - qir * bn[4]);
      real term6 = 4 * bn[3];


      real3 frc0;
      frc0.x = de * dr.x + term1 * di.x + term2 * dk.x +
         term3 * (diqkx - dkqix) + term4 * qix + term5 * qkx +
         term6 * (qixk + qkxi);
      frc0.y = de * dr.y + term1 * di.y + term2 * dk.y +
         term3 * (diqky - dkqiy) + term4 * qiy + term5 * qky +
         term6 * (qiyk + qkyi);
      frc0.z = de * dr.z + term1 * di.z + term2 * dk.z +
         term3 * (diqkz - dkqiz) + term4 * qiz + term5 * qkz +
         term6 * (qizk + qkzi);
      frc += (mscale * f) * frc0;

      // torque
      real dirx = di.y * dr.z - di.z * dr.y; // Di x r
      real diry = di.z * dr.x - di.x * dr.z;
      real dirz = di.x * dr.y - di.y * dr.x;
      real dkrx = dk.y * dr.z - dk.z * dr.y; // Dk x r
      real dkry = dk.z * dr.x - dk.x * dr.z;
      real dkrz = dk.x * dr.y - dk.y * dr.x;
      real dikx = di.y * dk.z - di.z * dk.y; // Di x Dk
      real diky = di.z * dk.x - di.x * dk.z;
      real dikz = di.x * dk.y - di.y * dk.x;

      real qirx = qiz * dr.y - qiy * dr.z; // r x (Qi r)
      real qiry = qix * dr.z - qiz * dr.x;
      real qirz = qiy * dr.x - qix * dr.y;
      real qkrx = qkz * dr.y - qky * dr.z; // r x (Qk r)
      real qkry = qkx * dr.z - qkz * dr.x;
      real qkrz = qky * dr.x - qkx * dr.y;
      real qikx = qky * qiz - qkz * qiy; // (Qk r) x (Qi r)
      real qiky = qkz * qix - qkx * qiz;
      real qikz = qkx * qiy - qky * qix;

      real qikrx = qizk * dr.y - qiyk * dr.z;
      real qikry = qixk * dr.z - qizk * dr.x;
      real qikrz = qiyk * dr.x - qixk * dr.y;
      real qkirx = qkzi * dr.y - qkyi * dr.z;
      real qkiry = qkxi * dr.z - qkzi * dr.x;
      real qkirz = qkyi * dr.x - qkxi * dr.y;

      real diqkrx = diqkz * dr.y - diqky * dr.z;
      real diqkry = diqkx * dr.z - diqkz * dr.x;
      real diqkrz = diqky * dr.x - diqkx * dr.y;
      real dkqirx = dkqiz * dr.y - dkqiy * dr.z;
      real dkqiry = dkqix * dr.z - dkqiz * dr.x;
      real dkqirz = dkqiy * dr.x - dkqix * dr.y;

      real dqikx = di.y * qkz - di.z * qky + dk.y * qiz - dk.z * qiy -
         2 *
            (qixy * qkxz + qiyy * qkyz + qiyz * qkzz - qixz * qkxy -
             qiyz * qkyy - qizz * qkyz);
      real dqiky = di.z * qkx - di.x * qkz + dk.z * qix - dk.x * qiz -
         2 *
            (qixz * qkxx + qiyz * qkxy + qizz * qkxz - qixx * qkxz -
             qixy * qkyz - qixz * qkzz);
      real dqikz = di.x * qky - di.y * qkx + dk.x * qiy - dk.y * qix -
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
      trq1 += (mscale * f) * trq0;
      trq0.x = bn[1] * dikx + term2 * dkrx - term3 * (dqikx + diqkrx) -
         term5 * qkrx - term6 * (qkirx - qikx);
      trq0.y = bn[1] * diky + term2 * dkry - term3 * (dqiky + diqkry) -
         term5 * qkry - term6 * (qkiry - qiky);
      trq0.z = bn[1] * dikz + term2 * dkrz - term3 * (dqikz + diqkrz) -
         term5 * qkrz - term6 * (qkirz - qikz);
      trq2 += (mscale * f) * trq0;
   }

   // epolar
   if CONSTEXPR (do_g) {
      f *= 0.5f;

      real uird = dot3(uid, dr);
      real ukrd = dot3(ukd, dr);
      real uirp = dot3(uip, dr);
      real ukrp = dot3(ukp, dr);

      real di_ukd = dot3(di, ukd);
      real di_ukp = dot3(di, ukp);
      real dk_uid = dot3(dk, uid);
      real dk_uip = dot3(dk, uip);

      real uid_qkr = dot3(uid, qk_dr); // <uid Qk r>
      real uip_qkr = dot3(uip, qk_dr); // <uip Qk r>
      real ukd_qir = dot3(ukd, qi_dr); // <ukd Qi r>
      real ukp_qir = dot3(ukp, qi_dr); // <ukp Qi r>

      // |Qi ukd>, |Qk uid>, |Qi ukp>, |Qk uip>
      real3 qi_ukd = matvec(qixx, qixy, qixz, qiyy, qiyz, qizz, ukd);
      real3 qk_uid = matvec(qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, uid);
      real3 qi_ukp = matvec(qixx, qixy, qixz, qiyy, qiyz, qizz, ukp);
      real3 qk_uip = matvec(qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, uip);


      real3 ufldi, ufldk;
      real dufldi[6], dufldk[6];

      // get the induced dipole field used for dipole torques
      ufldi = f * sr3 * (pscale * ukd + dscale * ukp) -
         f * sr5 * (pscale * ukrd + dscale * ukrp) * dr;
      ufldk = f * sr3 * (pscale * uid + dscale * uip) -
         f * sr5 * (pscale * uird + dscale * uirp) * dr;


      // get induced dipole field gradient used for quadrupole torques
      real tix5, tiy5, tiz5, tkx5, tky5, tkz5, tuir, tukr;
      tix5 = 2 * f * sr5 * (pscale * ukd.x + dscale * ukp.x);
      tiy5 = 2 * f * sr5 * (pscale * ukd.y + dscale * ukp.y);
      tiz5 = 2 * f * sr5 * (pscale * ukd.z + dscale * ukp.z);
      tkx5 = 2 * f * sr5 * (pscale * uid.x + dscale * uip.x);
      tky5 = 2 * f * sr5 * (pscale * uid.y + dscale * uip.y);
      tkz5 = 2 * f * sr5 * (pscale * uid.z + dscale * uip.z);
      tuir = -f * sr7 * (pscale * ukrd + dscale * ukrp);
      tukr = -f * sr7 * (pscale * uird + dscale * uirp);
      dufldi[0] = (dr.x * tix5 + dr.x * dr.x * tuir);
      dufldi[1] = (dr.x * tiy5 + dr.y * tix5 + 2 * dr.x * dr.y * tuir);
      dufldi[2] = (dr.y * tiy5 + dr.y * dr.y * tuir);
      dufldi[3] = (dr.x * tiz5 + dr.z * tix5 + 2 * dr.x * dr.z * tuir);
      dufldi[4] = (dr.y * tiz5 + dr.z * tiy5 + 2 * dr.y * dr.z * tuir);
      dufldi[5] = (dr.z * tiz5 + dr.z * dr.z * tuir);
      dufldk[0] = (-dr.x * tkx5 - dr.x * dr.x * tukr);
      dufldk[1] = (-dr.x * tky5 - dr.y * tkx5 - 2 * dr.x * dr.y * tukr);
      dufldk[2] = (-dr.y * tky5 - dr.y * dr.y * tukr);
      dufldk[3] = (-dr.x * tkz5 - dr.z * tkx5 - 2 * dr.x * dr.z * tukr);
      dufldk[4] = (-dr.y * tkz5 - dr.z * tky5 - 2 * dr.y * dr.z * tukr);
      dufldk[5] = (-dr.z * tkz5 - dr.z * dr.z * tukr);

      real3 trq0;
      trq0.x = di.z * ufldi.y - di.y * ufldi.z + qixz * dufldi[1] -
         qixy * dufldi[3] + 2 * qiyz * (dufldi[2] - dufldi[5]) +
         (qizz - qiyy) * dufldi[4];
      trq0.y = di.x * ufldi.z - di.z * ufldi.x - qiyz * dufldi[1] +
         qixy * dufldi[4] + 2 * qixz * (dufldi[5] - dufldi[0]) +
         (qixx - qizz) * dufldi[3];
      trq0.z = di.y * ufldi.x - di.x * ufldi.y + qiyz * dufldi[3] -
         qixz * dufldi[4] + 2 * qixy * (dufldi[0] - dufldi[2]) +
         (qiyy - qixx) * dufldi[1];
      trq1 += trq0;
      trq0.x = dk.z * ufldk.y - dk.y * ufldk.z + qkxz * dufldk[1] -
         qkxy * dufldk[3] + 2 * qkyz * (dufldk[2] - dufldk[5]) +
         (qkzz - qkyy) * dufldk[4];
      trq0.y = dk.x * ufldk.z - dk.z * ufldk.x - qkyz * dufldk[1] +
         qkxy * dufldk[4] + 2 * qkxz * (dufldk[5] - dufldk[0]) +
         (qkxx - qkzz) * dufldk[3];
      trq0.z = dk.y * ufldk.x - dk.x * ufldk.y + qkyz * dufldk[3] -
         qkxz * dufldk[4] + 2 * qkxy * (dufldk[0] - dufldk[2]) +
         (qkyy - qkxx) * dufldk[1];
      trq2 += trq0;

      // get the field gradient for direct polarization force
      real3 frcd = make_real3(0, 0, 0);
      real3 frcp = make_real3(0, 0, 0);
      // uind/p - charge
      frcd += sr3 * (ck * uip - ci * ukp) - sr5 * (ck * uirp - ci * ukrp) * dr;
      frcp += sr3 * (ck * uid - ci * ukd) - sr5 * (ck * uird - ci * ukrd) * dr;
      // uind/p - dipole
      frcd -= sr5 *
         (uirp * dk + ukrp * di + dir * ukp + dkr * uip +
          (di_ukp + dk_uip) * dr);
      frcp -= sr5 *
         (uird * dk + ukrd * di + dir * ukd + dkr * uid +
          (di_ukd + dk_uid) * dr);
      frcd += sr7 * (dir * ukrp + dkr * uirp) * dr;
      frcp += sr7 * (dir * ukrd + dkr * uird) * dr;
      // uind/p - quadrupole
      frcd +=
         2 * sr5 * (qi_ukp - qk_uip) + sr9 * (qir * ukrp - qkr * uirp) * dr;
      frcp +=
         2 * sr5 * (qi_ukd - qk_uid) + sr9 * (qir * ukrd - qkr * uird) * dr;
      frcd += 2 * sr7 * (uirp * qk_dr - ukrp * qi_dr) +
         2 * sr7 * (uip_qkr - ukp_qir) * dr + sr7 * (qkr * uip - qir * ukp);
      frcp += 2 * sr7 * (uird * qk_dr - ukrd * qi_dr) +
         2 * sr7 * (uid_qkr - ukd_qir) * dr + sr7 * (qkr * uid - qir * ukd);
      frc -= f * (dscale * frcd + pscale * frcp);

      // get the dtau/dr terms used for mutual polarization force
      real uid_ukp = dot3(uid, ukp);
      real uip_ukd = dot3(uip, ukd);
      frcd = sr5 * (uird * ukp + ukrd * uip + uirp * ukd + ukrp * uid);
      frcd += sr5 * (uid_ukp + uip_ukd) * dr;
      frcd -= sr7 * (uird * ukrp + ukrd * uirp) * dr;
      frc += f * uscale * frcd;
   }

   // save the results
   if CONSTEXPR (do_g) {
      frci += frc;
      frck -= frc;
      trqi += trq1;
      trqk += trq2;
   }
   if CONSTEXPR (do_v) {
      vtlxx -= dr.x * frc.x;
      vtlxy -= 0.5f * (dr.y * frc.x + dr.x * frc.y);
      vtlxz -= 0.5f * (dr.z * frc.x + dr.x * frc.z);
      vtlyy -= dr.y * frc.y;
      vtlyz -= 0.5f * (dr.z * frc.y + dr.y * frc.z);
      vtlzz -= dr.z * frc.z;
   }
}


// Rt Q = G
__device__
void rotQI2GVector(const real (&restrict rot)[3][3], real3 qif,
                   real3& restrict glf)
{
   glf = make_real3(dot3(rot[0][0], rot[1][0], rot[2][0], qif),
                    dot3(rot[0][1], rot[1][1], rot[2][1], qif),
                    dot3(rot[0][2], rot[1][2], rot[2][2], qif));
}


// R G = Q
__device__
void rotG2QIVector(const real (&restrict rot)[3][3], real3 glf,
                   real3& restrict qif)
{
   qif = make_real3(dot3(rot[0][0], rot[0][1], rot[0][2], glf),
                    dot3(rot[1][0], rot[1][1], rot[1][2], glf),
                    dot3(rot[2][0], rot[2][1], rot[2][2], glf));
}


// R G Rt = Q
__device__
void rotG2QIMat_v1(const real (&restrict rot)[3][3], //
                   real glxx, real glxy, real glxz,  //
                   real glyy, real glyz, real glzz,  //
                   real& restrict qixx, real& restrict qixy,
                   real& restrict qixz, real& restrict qiyy,
                   real& restrict qiyz, real& restrict qizz)
{
   real gl[3][3] = {{glxx, glxy, glxz}, {glxy, glyy, glyz}, {glxz, glyz, glzz}};
   real out[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
   // out[i][j] = sum(k,m) R[i][k] gl[k][m] Rt[m][j]
   //           = sum(k,m) R[i][k] gl[k][m] R[j][m]
   for (int i = 0; i < 2; ++i)
      for (int j = i; j < 3; ++j)
         for (int k = 0; k < 3; ++k)
            for (int m = 0; m < 3; ++m)
               out[i][j] += rot[i][k] * gl[k][m] * rot[j][m];
   qixx = out[0][0];
   qixy = out[0][1];
   qixz = out[0][2];
   qiyy = out[1][1];
   qiyz = out[1][2];
   // qizz = out[2][2];
   qizz = -(out[0][0] + out[1][1]);
}


// R G Rt = Q
__device__
void rotG2QIMat_v2(const real (&restrict r)[3][3],  //
                   real glxx, real glxy, real glxz, //
                   real glyy, real glyz, real glzz, //
                   real& restrict qixx, real& restrict qixy,
                   real& restrict qixz, real& restrict qiyy,
                   real& restrict qiyz, real& restrict qizz)
{
   // clang-format off
   qixx=r[0][0]*(r[0][0]*glxx+2*r[0][1]*glxy) + r[0][1]*(r[0][1]*glyy+2*r[0][2]*glyz) + r[0][2]*(r[0][2]*glzz+2*r[0][0]*glxz);
   qiyy=r[1][0]*(r[1][0]*glxx+2*r[1][1]*glxy) + r[1][1]*(r[1][1]*glyy+2*r[1][2]*glyz) + r[1][2]*(r[1][2]*glzz+2*r[1][0]*glxz);
   qixy=r[0][0]*(r[1][0]*glxx+r[1][1]*glxy+r[1][2]*glxz) + r[0][1]*(r[1][0]*glxy+r[1][1]*glyy+r[1][2]*glyz) + r[0][2]*(r[1][0]*glxz+r[1][1]*glyz+r[1][2]*glzz);
   qixz=r[0][0]*(r[2][0]*glxx+r[2][1]*glxy+r[2][2]*glxz) + r[0][1]*(r[2][0]*glxy+r[2][1]*glyy+r[2][2]*glyz) + r[0][2]*(r[2][0]*glxz+r[2][1]*glyz+r[2][2]*glzz);
   qiyz=r[1][0]*(r[2][0]*glxx+r[2][1]*glxy+r[2][2]*glxz) + r[1][1]*(r[2][0]*glxy+r[2][1]*glyy+r[2][2]*glyz) + r[1][2]*(r[2][0]*glxz+r[2][1]*glyz+r[2][2]*glzz);
   // clang-format on
   qizz = -(qixx + qiyy);
}


#define rotG2QIMatrix rotG2QIMat_v2


template <class Ver, class ETYP>
__device__
void pair_mplar_v2(                                                       //
   real r2, real3 dR, real mscale, real dscale, real pscale, real uscale, //
   real ci, real3 Id, real Iqxx, real Iqxy, real Iqxz, real Iqyy, real Iqyz,
   real Iqzz, real3 Iud, real3 Iup, real pdi, real pti, //
   real ck, real3 Kd, real Kqxx, real Kqxy, real Kqxz, real Kqyy, real Kqyz,
   real Kqzz, real3 Kud, real3 Kup, real pdk, real ptk, //
   real f, real aewald,                                 //
   real3& restrict frci, real3& restrict frck, real3& restrict trqi,
   real3& restrict trqk, real& restrict etl, real& restrict vtlxx,
   real& restrict vtlxy, real& restrict vtlxz, real& restrict vtlyy,
   real& restrict vtlyz, real& restrict vtlzz)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   // a rotation matrix that rotates (xr,yr,zr) to (0,0,r); R G = Q
   real rot[3][3];
   real bn[6];
   real sr3, sr5, sr7, sr9;
   real r = REAL_SQRT(r2);
   real invr1 = REAL_RECIP(r);
   {
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


      if CONSTEXPR (eq<ETYP, EWALD>()) {
         if CONSTEXPR (!do_g) {
            damp_ewald<5>(bn, r, invr1, rr2, aewald);
         } else {
            damp_ewald<6>(bn, r, invr1, rr2, aewald);
         }
      } else if CONSTEXPR (eq<ETYP, NON_EWALD>()) {
         bn[0] = rr1;
         bn[1] = rr3;
         bn[2] = rr5;
         bn[3] = rr7;
         bn[4] = rr9;
         if CONSTEXPR (do_g) {
            bn[5] = rr11;
         }
      }


      // if use_thole
      real ex3, ex5, ex7, ex9;
      damp_thole4(r, pdi, pti, pdk, ptk, ex3, ex5, ex7, ex9);
      sr3 = bn[1] - ex3 * rr3;
      sr5 = bn[2] - ex5 * rr5;
      sr7 = bn[3] - ex7 * rr7;
      sr9 = bn[4] - ex9 * rr9;
      // end if use_thole


      real3 rotz = invr1 * dR;
      // pick a random vector as rotx; rotx and rotz cannot be parallel
      real3 rotx = rotz;
      if (dR.y != 0 || dR.z != 0)
         rotx.x += 1;
      else
         rotx.y += 1;
      // Gram–Schmidt process for rotx with respect to rotz
      rotx -= dot3(rotx, rotz) * rotz;
      // normalize rotx
      real invxlen = REAL_RSQRT(dot3(rotx, rotx));
      rotx = invxlen * rotx;
      real3 roty = cross(rotz, rotx);
      rot[0][0] = rotx.x;
      rot[0][1] = rotx.y;
      rot[0][2] = rotx.z;
      rot[1][0] = roty.x;
      rot[1][1] = roty.y;
      rot[1][2] = roty.z;
      rot[2][0] = rotz.x;
      rot[2][1] = rotz.y;
      rot[2][2] = rotz.z;
   }


   real3 di, dk;
   rotG2QIVector(rot, Id, di);
   rotG2QIVector(rot, Kd, dk);
   real qixx, qixy, qixz, qiyy, qiyz, qizz;
   real qkxx, qkxy, qkxz, qkyy, qkyz, qkzz;
   rotG2QIMatrix(rot, Iqxx, Iqxy, Iqxz, Iqyy, Iqyz, Iqzz, qixx, qixy, qixz,
                 qiyy, qiyz, qizz);
   rotG2QIMatrix(rot, Kqxx, Kqxy, Kqxz, Kqyy, Kqyz, Kqzz, qkxx, qkxy, qkxz,
                 qkyy, qkyz, qkzz);
   real3 uid, uip;
   rotG2QIVector(rot, Iud, uid);
   rotG2QIVector(rot, Iup, uip);
   real3 ukd, ukp;
   rotG2QIVector(rot, Kud, ukd);
   rotG2QIVector(rot, Kup, ukp);


   // phi,dphi/d(x,y,z),d2phi/dd(xx,yy,zz,xy,xz,yz)
   //   0        1 2 3            4  5  6  7  8  9
   real phi1[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
   real phi2[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
   real phi1z[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};


   if CONSTEXPR (eq<ETYP, EWALD>()) {
      mscale = 1;
      dscale = 0.5f;
      pscale = 0.5f;
      uscale = 0.5f;
   } else {
      dscale *= 0.5f;
      pscale *= 0.5f;
      uscale *= 0.5f;
   }


   // C-C
   {
      real coef1 = bn[0];
      real coef3 = bn[1] * r;
      // phi_c c
      phi1[0] += coef1 * ck;
      phi2[0] += coef1 * ci;
      phi1z[0] += coef3 * ck;
   }


   // D-C and C-D
   {
      real coef3 = bn[1] * r;
      real coef5 = (bn[1] - bn[2] * r2);
      // phi_d c
      phi1[0] += -coef3 * dk.z;
      phi2[0] += coef3 * di.z;
      phi1z[0] += coef5 * dk.z;
      // dphi_c d
      // phi1[1]; phi1[2];
      phi1[3] += coef3 * ck;
      // phi2[1]; phi2[2];
      phi2[3] += -coef3 * ci;
      // phi1z[1]; phi1z[2];
      phi1z[3] += -coef5 * ck;
   }


   // D-D
   {
      real coef3 = bn[1];
      real coef5 = (bn[1] - bn[2] * r2);
      real coez5 = bn[2] * r;
      real coez7 = (3 * bn[2] - bn[3] * r2) * r;
      // dphi_d d
      phi1[1] += coef3 * dk.x;
      phi1[2] += coef3 * dk.y;
      phi1[3] += coef5 * dk.z;
      phi2[1] += coef3 * di.x;
      phi2[2] += coef3 * di.y;
      phi2[3] += coef5 * di.z;
      phi1z[1] += coez5 * dk.x;
      phi1z[2] += coez5 * dk.y;
      phi1z[3] += coez7 * dk.z;
   }


   // Q-C and C-Q
   {
      real coef3 = bn[1];
      real coef5 = bn[2] * r2;
      real coez5 = bn[2] * r;
      real coez7 = bn[3] * r2 * r;
      // phi_q c
      phi1[0] += coef5 * qkzz;
      phi2[0] += coef5 * qizz;
      phi1z[0] += -(2 * coez5 - coez7) * qkzz;
      // d2phi_c q
      phi1[4] += -coef3 * ck;
      phi1[5] += -coef3 * ck;
      phi1[6] += -(coef3 - coef5) * ck;
      // phi1[7]; phi1[8]; phi1[9];
      phi2[4] += -coef3 * ci;
      phi2[5] += -coef3 * ci;
      phi2[6] += -(coef3 - coef5) * ci;
      // phi2[7]; phi2[8]; phi2[9];
      phi1z[4] += -coez5 * ck;
      phi1z[5] += -coez5 * ck;
      phi1z[6] += -(3 * coez5 - coez7) * ck;
      // phi1z[7]; phi1z[8]; phi1z[9];
   }


   // Q-D and D-Q
   {
      real coef5 = bn[2] * r;
      real coef7 = bn[3] * r2 * r;
      real coez7 = (bn[2] - bn[3] * r2);
      real coez9 = (3 * bn[3] - bn[4] * r2) * r2;
      // dphi_q d
      phi1[1] += -2 * coef5 * qkxz;
      phi1[2] += -2 * coef5 * qkyz;
      phi1[3] += -(2 * coef5 - coef7) * qkzz;
      phi2[1] += 2 * coef5 * qixz;
      phi2[2] += 2 * coef5 * qiyz;
      phi2[3] += (2 * coef5 - coef7) * qizz;
      phi1z[1] += 2 * coez7 * qkxz;
      phi1z[2] += 2 * coez7 * qkyz;
      phi1z[3] += (2 * coez7 - coez9) * qkzz;
      // d2phi_d q
      phi1[4] += coef5 * dk.z;
      phi1[5] += coef5 * dk.z;
      phi1[6] += (3 * coef5 - coef7) * dk.z;
      // phi1[7];
      phi1[8] += 2 * coef5 * dk.x;
      phi1[9] += 2 * coef5 * dk.y;
      //
      phi2[4] += -coef5 * di.z;
      phi2[5] += -coef5 * di.z;
      phi2[6] += -(3 * coef5 - coef7) * di.z;
      // phi2[7];
      phi2[8] += -2 * coef5 * di.x;
      phi2[9] += -2 * coef5 * di.y;
      //
      phi1z[4] += -coez7 * dk.z;
      phi1z[5] += -coez7 * dk.z;
      phi1z[6] += -(3 * coez7 - coez9) * dk.z;
      // phi1z[7];
      phi1z[8] += -2 * coez7 * dk.x;
      phi1z[9] += -2 * coez7 * dk.y;
   }


   // Q-Q
   {
      // d2phi_q q
      real coef5 = bn[2];
      real coef7 = bn[3] * r2;
      real coef9 = bn[4] * r2 * r2;
      real coez7 = bn[3] * r;
      real coez9 = bn[4] * r2 * r;
      real coez11 = bn[5] * r2 * r2 * r;
      //
      phi1[4] += 2 * coef5 * qkxx - coef7 * qkzz;
      phi1[5] += 2 * coef5 * qkyy - coef7 * qkzz;
      phi1[6] += (2 * coef5 - 5 * coef7 + coef9) * qkzz;
      phi1[7] += 4 * coef5 * qkxy;
      phi1[8] += 4 * (coef5 - coef7) * qkxz;
      phi1[9] += 4 * (coef5 - coef7) * qkyz;
      //
      phi2[4] += 2 * coef5 * qixx - coef7 * qizz;
      phi2[5] += 2 * coef5 * qiyy - coef7 * qizz;
      phi2[6] += (2 * coef5 - 5 * coef7 + coef9) * qizz;
      phi2[7] += 4 * coef5 * qixy;
      phi2[8] += 4 * (coef5 - coef7) * qixz;
      phi2[9] += 4 * (coef5 - coef7) * qiyz;
      //
      phi1z[4] += 2 * coez7 * qkxx + (2 * coez7 - coez9) * qkzz;
      phi1z[5] += 2 * coez7 * qkyy + (2 * coez7 - coez9) * qkzz;
      phi1z[6] += (12 * coez7 - 9 * coez9 + coez11) * qkzz;
      phi1z[7] += 4 * coez7 * qkxy;
      phi1z[8] += 4 * (3 * coez7 - coez9) * qkxz;
      phi1z[9] += 4 * (3 * coez7 - coez9) * qkyz;
   }


   #pragma unroll
   for (int i = 0; i < 10; ++i) {
      phi1[i] *= mscale;
      phi2[i] *= mscale;
      phi1z[i] *= mscale;
   }


   if CONSTEXPR (do_e) {
      real e = phi1[0] * ci + phi1[1] * di.x + phi1[2] * di.y + phi1[3] * di.z +
         phi1[4] * qixx + phi1[5] * qiyy + phi1[6] * qizz + phi1[7] * qixy +
         phi1[8] * qixz + phi1[9] * qiyz;
      etl += f * e;
   }


   real phi1d[3] = {0, 0, 0};
   real phi2d[3] = {0, 0, 0};
   real phi1dz[3] = {0, 0, 0};


   // U-C and C-U
   {
      real coe3 = sr3 * r;
      real coe5 = sr3 - sr5 * r2;
      real coed3 = dscale * coe3;
      real coed5 = dscale * coe5;
      real coep3 = pscale * coe3;
      real coep5 = pscale * coe5;
      // phi_u c
      phi1[0] += -(coed3 * ukp.z + coep3 * ukd.z);
      phi2[0] += coed3 * uip.z + coep3 * uid.z;
      phi1z[0] += coed5 * ukp.z + coep5 * ukd.z;
      // dphi_c u
      phi1d[2] += coe3 * ck;
      phi2d[2] += -coe3 * ci;
      phi1dz[2] += -coe5 * ck;
   }


   // U-D and D-U
   {
      real coe3 = sr3;
      real coe5 = sr5 * r2;
      real coez5 = sr5 * r;
      real coez7 = sr7 * r2 * r;
      real coed3 = dscale * coe3;
      real coed5 = dscale * coe5;
      real coedz5 = dscale * coez5;
      real coedz7 = dscale * coez7;
      real coep3 = pscale * coe3;
      real coep5 = pscale * coe5;
      real coepz5 = pscale * coez5;
      real coepz7 = pscale * coez7;
      // dphi_u d
      phi1[1] += coed3 * ukp.x + coep3 * ukd.x;
      phi1[2] += coed3 * ukp.y + coep3 * ukd.y;
      phi1[3] += (coed3 - coed5) * ukp.z + (coep3 - coep5) * ukd.z;
      phi2[1] += coed3 * uip.x + coep3 * uid.x;
      phi2[2] += coed3 * uip.y + coep3 * uid.y;
      phi2[3] += (coed3 - coed5) * uip.z + (coep3 - coep5) * uid.z;
      phi1z[1] += coedz5 * ukp.x + coepz5 * ukd.x;
      phi1z[2] += coedz5 * ukp.y + coepz5 * ukd.y;
      phi1z[3] += (3 * coedz5 - coedz7) * ukp.z + (3 * coepz5 - coepz7) * ukd.z;
      // dphi_d u
      phi1d[0] += coe3 * dk.x;
      phi1d[1] += coe3 * dk.y;
      phi1d[2] += (coe3 - coe5) * dk.z;
      phi2d[0] += coe3 * di.x;
      phi2d[1] += coe3 * di.y;
      phi2d[2] += (coe3 - coe5) * di.z;
      phi1dz[0] += coez5 * dk.x;
      phi1dz[1] += coez5 * dk.y;
      phi1dz[2] += (3 * coez5 - coez7) * dk.z;
   }


   // U-Q and Q-U
   {
      real coe5 = sr5 * r;
      real coe7 = sr7 * r2 * r;
      real coez7 = sr5 - sr7 * r2;
      real coez9 = (3 * sr7 - sr9 * r2) * r2;
      real coed5 = dscale * coe5;
      real coed7 = dscale * coe7;
      real coedz7 = dscale * coez7;
      real coedz9 = dscale * coez9;
      real coep5 = pscale * coe5;
      real coep7 = pscale * coe7;
      real coepz7 = pscale * coez7;
      real coepz9 = pscale * coez9;
      // d2phi_u q
      phi1[4] += coed5 * ukp.z + coep5 * ukd.z;
      phi1[5] += coed5 * ukp.z + coep5 * ukd.z;
      phi1[6] += (3 * coed5 - coed7) * ukp.z + (3 * coep5 - coep7) * ukd.z;
      // phi1[7];
      phi1[8] += 2 * (coed5 * ukp.x + coep5 * ukd.x);
      phi1[9] += 2 * (coed5 * ukp.y + coep5 * ukd.y);
      //
      phi2[4] += -(coed5 * uip.z + coep5 * uid.z);
      phi2[5] += -(coed5 * uip.z + coep5 * uid.z);
      phi2[6] += -(3 * coed5 - coed7) * uip.z - (3 * coep5 - coep7) * uid.z;
      // phi2[7];
      phi2[8] += -2 * (coed5 * uip.x + coep5 * uid.x);
      phi2[9] += -2 * (coed5 * uip.y + coep5 * uid.y);
      //
      phi1z[4] += -(coedz7 * ukp.z + coepz7 * ukd.z);
      phi1z[5] += -(coedz7 * ukp.z + coepz7 * ukd.z);
      phi1z[6] +=
         -(3 * coedz7 - coedz9) * ukp.z - (3 * coepz7 - coepz9) * ukd.z;
      // phi1z[7];
      phi1z[8] += -2 * (coedz7 * ukp.x + coepz7 * ukd.x);
      phi1z[9] += -2 * (coedz7 * ukp.y + coepz7 * ukd.y);
      // dphi_q u
      phi1d[0] += -2 * coe5 * qkxz;
      phi1d[1] += -2 * coe5 * qkyz;
      phi1d[2] += -(2 * coe5 - coe7) * qkzz;
      phi2d[0] += 2 * coe5 * qixz;
      phi2d[1] += 2 * coe5 * qiyz;
      phi2d[2] += (2 * coe5 - coe7) * qizz;
      phi1dz[0] += 2 * coez7 * qkxz;
      phi1dz[1] += 2 * coez7 * qkyz;
      phi1dz[2] += (2 * coez7 - coez9) * qkzz;
   }


   real3 frc, trq1, trq2;
   if CONSTEXPR (do_g) {
      // torque
      real3 trqa = cross(phi1[1], phi1[2], phi1[3], di);
      trqa.x += phi1[9] * (qizz - qiyy) + 2 * (phi1[5] - phi1[6]) * qiyz +
         phi1[7] * qixz - phi1[8] * qixy;
      trqa.y += phi1[8] * (qixx - qizz) + 2 * (phi1[6] - phi1[4]) * qixz +
         phi1[9] * qixy - phi1[7] * qiyz;
      trqa.z += phi1[7] * (qiyy - qixx) + 2 * (phi1[4] - phi1[5]) * qixy +
         phi1[8] * qiyz - phi1[9] * qixz;
      real3 trqb = cross(phi2[1], phi2[2], phi2[3], dk);
      trqb.x += phi2[9] * (qkzz - qkyy) + 2 * (phi2[5] - phi2[6]) * qkyz +
         phi2[7] * qkxz - phi2[8] * qkxy;
      trqb.y += phi2[8] * (qkxx - qkzz) + 2 * (phi2[6] - phi2[4]) * qkxz +
         phi2[9] * qkxy - phi2[7] * qkyz;
      trqb.z += phi2[7] * (qkyy - qkxx) + 2 * (phi2[4] - phi2[5]) * qkxy +
         phi2[8] * qkyz - phi2[9] * qkxz;
      trq1 = trqa;
      trq2 = trqb;


      real3 trqau =
         cross(phi1d[0], phi1d[1], phi1d[2], (dscale * uip + pscale * uid));
      real3 trqbu =
         cross(phi2d[0], phi2d[1], phi2d[2], (dscale * ukp + pscale * ukd));


      // gradient
      real frc1z = phi1z[0] * ci + phi1z[1] * di.x + phi1z[2] * di.y +
         phi1z[3] * di.z + phi1z[4] * qixx + phi1z[5] * qiyy + phi1z[6] * qizz +
         phi1z[7] * qixy + phi1z[8] * qixz + phi1z[9] * qiyz;
      frc1z +=
         dot3(phi1dz[0], phi1dz[1], phi1dz[2], (dscale * uip + pscale * uid));
      frc.x = -invr1 * (trqa.y + trqb.y + trqau.y + trqbu.y);
      frc.y = invr1 * (trqa.x + trqb.x + trqau.x + trqbu.x);
      frc.z = frc1z;
   }


   // U-U
   {
      real coeu5 = uscale * sr5 * r;
      real coeu7 = uscale * sr7 * r2 * r;
      frc.x += coeu5 *
         (uid.x * ukp.z + uid.z * ukp.x + uip.x * ukd.z + uip.z * ukd.x);
      frc.y += coeu5 *
         (uid.y * ukp.z + uid.z * ukp.y + uip.y * ukd.z + uip.z * ukd.y);
      frc.z += coeu5 *
            (uid.x * ukp.x + uid.y * ukp.y + uip.x * ukd.x + uip.y * ukd.y) +
         (3 * coeu5 - coeu7) * (uid.z * ukp.z + uip.z * ukd.z);
   }


   if CONSTEXPR (do_g) {
      real3 glfrc;
      rotQI2GVector(rot, frc, glfrc);
      frc = f * glfrc;
      frci += frc;
      frck -= frc;
      real3 gltrq1;
      rotQI2GVector(rot, trq1, gltrq1);
      trqi += f * gltrq1;
      real3 gltrq2;
      rotQI2GVector(rot, trq2, gltrq2);
      trqk += f * gltrq2;
   }
   if CONSTEXPR (do_v) {
      vtlxx -= dR.x * frc.x;
      vtlxy -= 0.5f * (dR.y * frc.x + dR.x * frc.y);
      vtlxz -= 0.5f * (dR.z * frc.x + dR.x * frc.z);
      vtlyy -= dR.y * frc.y;
      vtlyz -= 0.5f * (dR.z * frc.y + dR.y * frc.z);
      vtlzz -= dR.z * frc.z;
   }
}


#define pair_mplar pair_mplar_v2


#define EMPLAR_ARGS                                                            \
   size_t bufsize, energy_buffer restrict ebuf,                                \
      virial_buffer restrict vir_ebuf, grad_prec *restrict gx,                 \
      grad_prec *restrict gy, grad_prec *restrict gz, real *restrict trqx,     \
      real *restrict trqy, real *restrict trqz, TINKER_IMAGE_PARAMS, real off, \
      real f, const real(*restrict rpole)[10], const real *restrict pdamp,     \
      const real *restrict thole, const real(*restrict uind)[3],               \
      const real(*restrict uinp)[3], real(*restrict ufld)[3],                  \
      real(*restrict dufld)[6]


template <class Ver, class ETYP>
__global__
void emplar_cu1(EMPLAR_ARGS, int n, const Spatial::SortedAtom* restrict sorted,
                int niak, const int* restrict iak, const int* restrict lst,
                real aewald)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;
   static_assert(!Ver::a, "");


   const int iwarp = (threadIdx.x + blockIdx.x * blockDim.x) / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);
   const int offset = (threadIdx.x + blockIdx.x * blockDim.x) & (bufsize - 1);


   struct Data
   {
      real3 pos;
      real3 frc, trq; // force and torque
      real c;         // charge
      real3 d;        // dipole
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
            if CONSTEXPR (eq<ETYP, EWALD>()) {
               pair_mplar<Ver, EWALD>(
                  r2, dr, 1, 1, 1, 1, //
                  idat.c, idat.d, idat.qxx, idat.qxy, idat.qxz, idat.qyy,
                  idat.qyz, idat.qzz, idat.ud, idat.up, idat.damp,
                  idat.thole, //
                  data[klane].c, data[klane].d, data[klane].qxx,
                  data[klane].qxy, data[klane].qxz, data[klane].qyy,
                  data[klane].qyz, data[klane].qzz, data[klane].ud,
                  data[klane].up, data[klane].damp, data[klane].thole, //
                  f, aewald,                                           //
                  idat.frc, data[klane].frc, idat.trq, data[klane].trq, etl,
                  vtlxx, vtlxy, vtlxz, vtlyy, vtlyz, vtlzz);
            }
            if CONSTEXPR (eq<ETYP, NON_EWALD>()) {
               pair_mplar<Ver, NON_EWALD>(
                  r2, dr, 1, 1, 1, 1, //
                  idat.c, idat.d, idat.qxx, idat.qxy, idat.qxz, idat.qyy,
                  idat.qyz, idat.qzz, idat.ud, idat.up, idat.damp,
                  idat.thole, //
                  data[klane].c, data[klane].d, data[klane].qxx,
                  data[klane].qxy, data[klane].qxz, data[klane].qyy,
                  data[klane].qyz, data[klane].qzz, data[klane].ud,
                  data[klane].up, data[klane].damp,
                  data[klane].thole, //
                  f, 0,              //
                  idat.frc, data[klane].frc, idat.trq, data[klane].trq, etl,
                  vtlxx, vtlxy, vtlxz, vtlyy, vtlyz, vtlzz);
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
      }
      if CONSTEXPR (do_v)
         atomic_add(vtlxx, vtlxy, vtlxz, vtlyy, vtlyz, vtlzz, vir_ebuf, offset);
   } // end for (iw)
}


template <class Ver>
__global__
void emplar_cu2(EMPLAR_ARGS, const real* restrict x, const real* restrict y,
                const real* restrict z, int nmdpuexclude,
                const int (*restrict mdpuexclude)[2],
                const real (*restrict mdpuexclude_scale)[4])
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;
   static_assert(!Ver::a, "");


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
         real3 frci, frck, trqi, trqk;
         if CONSTEXPR (do_g) {
            frci = make_real3(0, 0, 0);
            frck = make_real3(0, 0, 0);
            trqi = make_real3(0, 0, 0);
            trqk = make_real3(0, 0, 0);
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


         pair_mplar<Ver, NON_EWALD>(
            r2, make_real3(xr, yr, zr), mscale, dscale, pscale, uscale, //
            ci, make_real3(dix, diy, diz), qixx, qixy, qixz, qiyy, qiyz, qizz,
            make_real3(uix, uiy, uiz), make_real3(uixp, uiyp, uizp), pdi,
            pti, //
            ck, make_real3(dkx, dky, dkz), qkxx, qkxy, qkxz, qkyy, qkyz, qkzz,
            make_real3(ukx, uky, ukz), make_real3(ukxp, ukyp, ukzp), pdk,
            ptk,  //
            f, 0, //
            frci, frck, trqi, trqk, e, vxx, vxy, vxz, vyy, vyz, vzz);


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
         }
         if CONSTEXPR (do_v)
            atomic_add(vxx, vxy, vxz, vyy, vyz, vzz, vir_ebuf, offset);
      } // end if (include)
   }
}


template <class Ver, class ETYP>
void emplar_cu(const real (*uind)[3], const real (*uinp)[3])
{
   const auto& st = *mspatial_unit;
   real off;
   if CONSTEXPR (eq<ETYP, EWALD>())
      off = switch_off(switch_ewald);
   else
      off = switch_off(switch_mpole);
   auto bufsize = buffer_size();


   const real f = electric / dielec;
   real aewald = 0;
   if CONSTEXPR (eq<ETYP, EWALD>()) {
      assert(epme_unit == ppme_unit);
      PMEUnit pu = epme_unit;
      aewald = pu->aewald;


      if CONSTEXPR (Ver::e) {
         auto ker0 = empole_self_cu<Ver::a>;
         launch_k1s(nonblk, n, ker0, //
                    bufsize, nullptr, em, rpole, n, f, aewald);
      }
   }
   if (st.niak > 0) {
      auto ker1 = emplar_cu1<Ver, ETYP>;
      launch_k1s(nonblk, WARP_SIZE * st.niak, ker1, //
                 bufsize, em, vir_em, gx, gy, gz, trqx, trqy, trqz,
                 TINKER_IMAGE_ARGS, off, f, rpole, pdamp, thole, uind, uinp,
                 ufld, dufld, //
                 n, st.sorted, st.niak, st.iak, st.lst, aewald);
   }
   if (nmdpuexclude > 0) {
      auto ker2 = emplar_cu2<Ver>;
      launch_k1s(nonblk, nmdpuexclude, ker2, //
                 bufsize, em, vir_em, gx, gy, gz, trqx, trqy, trqz,
                 TINKER_IMAGE_ARGS, off, f, rpole, pdamp, thole, uind, uinp,
                 ufld, dufld, //
                 x, y, z, nmdpuexclude, mdpuexclude, mdpuexclude_scale);
   }
}


template <class Ver>
void emplar_ewald_cu()
{
   // induce
   induce(uind, uinp);

   // empole real self; epolar real without epolar energy
   emplar_cu<Ver, EWALD>(uind, uinp);
   // empole recip
   empole_ewald_recip(Ver::value);
   // epolar recip self; must toggle off the calc::energy flag
   epolar_ewald_recip_self(Ver::value & ~calc::energy);

   // epolar energy
   if CONSTEXPR (Ver::e)
      epolar0_dotprod(uind, udirp);
}


template <class Ver>
void emplar_nonewald_cu()
{
   // induce
   induce(uind, uinp);

   // empole and epolar
   emplar_cu<Ver, NON_EWALD>(uind, uinp);
   if CONSTEXPR (Ver::e)
      epolar0_dotprod(uind, udirp);
}


void emplar_cu(int vers)
{
   if (use_ewald()) {
      if (vers == calc::v0)
         emplar_ewald_cu<calc::V0>();
      else if (vers == calc::v1)
         emplar_ewald_cu<calc::V1>();
      // else if (vers == calc::v3)
      //    emplar_ewald_cu<calc::V3>();
      else if (vers == calc::v4)
         emplar_ewald_cu<calc::V4>();
      else if (vers == calc::v5)
         emplar_ewald_cu<calc::V5>();
      else if (vers == calc::v6)
         emplar_ewald_cu<calc::V6>();
   } else {
      if (vers == calc::v0)
         emplar_nonewald_cu<calc::V0>();
      else if (vers == calc::v1)
         emplar_nonewald_cu<calc::V1>();
      // else if (vers == calc::v3)
      //    emplar_nonewald_cu<calc::V3>();
      else if (vers == calc::v4)
         emplar_nonewald_cu<calc::V4>();
      else if (vers == calc::v5)
         emplar_nonewald_cu<calc::V5>();
      else if (vers == calc::v6)
         emplar_nonewald_cu<calc::V6>();
   }
}
TINKER_NAMESPACE_END
