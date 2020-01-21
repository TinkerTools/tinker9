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
void pair_mplar(                                                          //
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
   constexpr int do_e = USE & calc::energy;
   constexpr int do_g = USE & calc::grad;
   constexpr int do_v = USE & calc::virial;

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

      // if use_thole
      real ex3, ex5, ex7, ex9;
      damp_thole4(r, rr2, dr.x, dr.y, dr.z, pdi, pti, pdk, ptk, ex3, ex5, ex7,
                  ex9);
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
      // clang-format off
      ufldi = f*sr3*(pscale*ukd +dscale*ukp) -f*sr5*(pscale*ukrd +dscale*ukrp)*dr;
      ufldk = f*sr3*(pscale*uid +dscale*uip) -f*sr5*(pscale*uird +dscale*uirp)*dr;
      // clang-format on

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
      // clang-format off
      // uind/p - charge
      frcd += sr3*(ck*uip -ci*ukp) -sr5*(ck*uirp -ci*ukrp)*dr;
      frcp += sr3*(ck*uid -ci*ukd) -sr5*(ck*uird -ci*ukrd)*dr;
      // uind/p - dipole
      frcd -= sr5*(uirp*dk +ukrp*di +dir*ukp +dkr*uip +(di_ukp+dk_uip)*dr);
      frcp -= sr5*(uird*dk +ukrd*di +dir*ukd +dkr*uid +(di_ukd+dk_uid)*dr);
      frcd += sr7*(dir*ukrp +dkr*uirp)*dr;
      frcp += sr7*(dir*ukrd +dkr*uird)*dr;
      // uind/p - quadrupole
      frcd += 2*sr5*(qi_ukp -qk_uip) + sr9*(qir*ukrp -qkr*uirp)*dr;
      frcp += 2*sr5*(qi_ukd -qk_uid) + sr9*(qir*ukrd -qkr*uird)*dr;
      frcd += 2*sr7*(uirp*qk_dr -ukrp*qi_dr) +2*sr7*(uip_qkr -ukp_qir)*dr +sr7*(qkr*uip -qir*ukp);
      frcp += 2*sr7*(uird*qk_dr -ukrd*qi_dr) +2*sr7*(uid_qkr -ukd_qir)*dr +sr7*(qkr*uid -qir*ukd);
      // clang-format on
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
__global__
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
            if CONSTEXPR (ETYP == elec_t::ewald) {
               pair_mplar<USE, elec_t::ewald>(
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
            if CONSTEXPR (ETYP == elec_t::coulomb) {
               pair_mplar<USE, elec_t::coulomb>(
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


         pair_mplar<USE, elec_t::coulomb>(
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


   // empole real self
   // epolar real gradient
   emplar_tmpl_cu<USE, elec_t::ewald>(uind, uinp);
   // epolar torque
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


   // empole and epolar
   emplar_tmpl_cu<USE, elec_t::coulomb>(uind, uinp);
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
