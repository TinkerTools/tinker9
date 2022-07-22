#pragma once
#include "math/switch.h"
#include "seq/damp_hippo.h"
#include "ff/hippomod.h"

namespace tinker {
#pragma acc routine seq
SEQ_CUDA
inline void pair_alterpol(ExpolScr scrtyp, real r, real r2, real pscale, real cut, real off, real xr, real yr, real zr,
   real springi, real sizi, real alphai, real springk, real sizk, real alphak,
   real ks2i[3][3], real ks2k[3][3])
{
   real cut2 = cut * cut;
   real sizik = sizi * sizk;
   real s2;
   real ds2;

   damp_expl(scrtyp, s2, ds2, r, sizik, alphai, alphak);

   if (r2 > cut2) {
      real taper, dtaper;
      switchTaper5<0>(r, cut, off, taper, dtaper);
      s2 = s2 * taper;
   }
   real p33i, p33k;
   p33i = springi * s2 * pscale;
   p33k = springk * s2 * pscale;

   real ai[3], ak[3];

   ai[0] = xr / r;
   ai[1] = yr / r;
   ai[2] = zr / r;

   ak[0] = -ai[0];
   ak[1] = -ai[1];
   ak[2] = -ai[2];
   #pragma acc loop seq
   for (int i{0}; i < 3; ++i) {
      #pragma acc loop seq
      for (int j{0}; j < 3; ++j) {
         ks2i[j][i] = p33i * ai[i] * ai[j];
         ks2k[j][i] = p33k * ak[i] * ak[j];
      }
   }
}
}
