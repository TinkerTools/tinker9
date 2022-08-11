#pragma once
#include "ff/hippo/expolscr.h"
#include "math/switch.h"
#include "seq/seq.h"

namespace tinker {
#pragma acc routine seq
template <bool DO_G>
SEQ_CUDA
inline void damp_expl(ExpolScr scrtyp, real& restrict s2, real& restrict ds2, real r, real sizik,
   real alphai, real alphak)
{
   constexpr real inv2 = 1. / 2, inv3 = 1. / 3;
   constexpr real one = 1.;

   real alphaik, dmpik2, dampik, dampik2, expik, s;

   if (scrtyp == ExpolScr::S2U) {
      alphaik = REAL_SQRT(alphai * alphak);
      dmpik2 = inv2 * alphaik;
      dampik = dmpik2 * r;
      dampik2 = dampik * dampik;
      expik = REAL_EXP(-dampik);
      s = (one + dampik + dampik2 * inv3) * expik;
      s2 = s * s;
      if (DO_G)
         ds2 = s * (-alphaik * inv3) * (dampik + dampik2) * expik;
   }

   s2 = sizik * s2;
   if (DO_G)
      ds2 = sizik * ds2;
}

SEQ_ROUTINE
inline void pair_alterpol(ExpolScr scrtyp, real r, real pscale, real cut, real off, real xr,
   real yr, real zr, real springi, real sizi, real alphai, real springk, real sizk, real alphak,
   real ks2i[3][3], real ks2k[3][3])
{
   real sizik = sizi * sizk;
   real s2;
   real ds2;

   constexpr bool DO_G = false;
   damp_expl<DO_G>(scrtyp, s2, ds2, r, sizik, alphai, alphak);

   if (r > cut) {
      real taper, dtaper;
      switchTaper5<0>(r, cut, off, taper, dtaper);
      s2 = s2 * taper;
   }

   real p33i, p33k;
   p33i = springi * s2 * pscale;
   p33k = springk * s2 * pscale;

   real ai[3]; // ak = -ai

   ai[0] = xr / r;
   ai[1] = yr / r;
   ai[2] = zr / r;

   #pragma acc loop seq
   for (int i = 0; i < 3; ++i) {
      #pragma acc loop seq
      for (int j = 0; j < 3; ++j) {
         ks2i[j][i] = p33i * ai[i] * ai[j];
         ks2k[j][i] = p33k * ai[i] * ai[j]; // ak_i * ak_j = ai_i * ai_j
      }
   }
}

SEQ_ROUTINE
inline void pair_dexpol(ExpolScr scrtyp, real r, real pscale, real cut, real off, real xr, real yr,
   real zr, real uix, real uiy, real uiz, real ukx, real uky, real ukz, real springi, real sizi,
   real alphai, real springk, real sizk, real alphak, const real f, real frc[3])
{
   real sizik = sizi * sizk;
   real s2;
   real ds2;

   constexpr bool DO_G = true;
   damp_expl<DO_G>(scrtyp, s2, ds2, r, sizik, alphai, alphak);

   if (r > cut) {
      real taper, dtaper;
      switchTaper5<1>(r, cut, off, taper, dtaper);
      ds2 = ds2 * taper + s2 * dtaper;
      s2 = s2 * taper;
   }
   real s2i = springi * s2 * pscale;
   real s2k = springk * s2 * pscale;
   real ds2i = springi * ds2 * pscale;
   real ds2k = springk * ds2 * pscale;

   // compute rotation matrix
   real ai[3][3], ak[3][3];
   ai[0][2] = xr / r;
   ai[1][2] = yr / r;
   ai[2][2] = zr / r;
   xr = 1.;
   yr = 0.;
   zr = 0.;
   real dot = ai[0][2];
   real eps = 1. / REAL_SQRT(2);
   if (fabs(dot) > eps) {
      xr = 0.;
      yr = 1.;
      dot = ai[1][2];
   }
   xr = xr - dot * ai[0][2];
   yr = yr - dot * ai[1][2];
   zr = zr - dot * ai[2][2];
   real dr = REAL_SQRT(xr * xr + yr * yr + zr * zr);
   ai[0][0] = xr / dr;
   ai[1][0] = yr / dr;
   ai[2][0] = zr / dr;
   ai[0][1] = ai[2][0] * ai[1][2] - ai[1][0] * ai[2][2];
   ai[1][1] = ai[0][0] * ai[2][2] - ai[2][0] * ai[0][2];
   ai[2][1] = ai[1][0] * ai[0][2] - ai[0][0] * ai[1][2];
   ak[0][0] = ai[0][0];
   ak[1][0] = ai[1][0];
   ak[2][0] = ai[2][0];
   ak[0][1] = -ai[0][1];
   ak[1][1] = -ai[1][1];
   ak[2][1] = -ai[2][1];
   ak[0][2] = -ai[0][2];
   ak[1][2] = -ai[1][2];
   ak[2][2] = -ai[2][2];

   // local frame force
   real frcil[3], frckl[3];
   real uixl = uix * ai[0][0] + uiy * ai[1][0] + uiz * ai[2][0];
   real uiyl = uix * ai[0][1] + uiy * ai[1][1] + uiz * ai[2][1];
   real uizl = uix * ai[0][2] + uiy * ai[1][2] + uiz * ai[2][2];
   real ukxl = -(ukx * ak[0][0] + uky * ak[1][0] + ukz * ak[2][0]);
   real ukyl = -(ukx * ak[0][1] + uky * ak[1][1] + ukz * ak[2][1]);
   real ukzl = -(ukx * ak[0][2] + uky * ak[1][2] + ukz * ak[2][2]);
   frcil[2] = REAL_POW(uizl, 2) * ds2i;
   frckl[2] = REAL_POW(ukzl, 2) * ds2k;
   // local frame torque
   constexpr real two = 2.;
   real tqxil = two * uiyl * uizl * s2i;
   real tqyil = -two * uixl * uizl * s2i;
   real tqxkl = two * ukyl * ukzl * s2k;
   real tqykl = -two * ukxl * ukzl * s2k;
   // convert local frame torques to local frame forces
   frcil[0] = -tqyil / r;
   frcil[1] = tqxil / r;
   frckl[0] = -tqykl / r;
   frckl[1] = tqxkl / r;
   // rotate force to global frame
   real frcxi = ai[0][0] * frcil[0] + ai[0][1] * frcil[1] + ai[0][2] * frcil[2];
   real frcyi = ai[1][0] * frcil[0] + ai[1][1] * frcil[1] + ai[1][2] * frcil[2];
   real frczi = ai[2][0] * frcil[0] + ai[2][1] * frcil[1] + ai[2][2] * frcil[2];
   real frcxk = ak[0][0] * frckl[0] + ak[0][1] * frckl[1] + ak[0][2] * frckl[2];
   real frcyk = ak[1][0] * frckl[0] + ak[1][1] * frckl[1] + ak[1][2] * frckl[2];
   real frczk = ak[2][0] * frckl[0] + ak[2][1] * frckl[1] + ak[2][2] * frckl[2];
   frc[0] = f * (frcxk - frcxi);
   frc[1] = f * (frcyk - frcyi);
   frc[2] = f * (frczk - frczi);
}
}
