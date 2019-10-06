#pragma once

#include "mathfunc.h"

TINKER_NAMESPACE_BEGIN
#pragma acc routine seq
CUDA_DEVICE_FUNCTION
inline void damp_thole2 (real r, real pdi, real pti, real pdk, real ptk,
                         real& RESTRICT scale3, real& RESTRICT scale5)
{
   scale3 = 1;
   scale5 = 1;
   real damp = pdi * pdk;
   if (damp != 0) {
      real pgamma = REAL_MIN (pti, ptk);
      damp = -pgamma * REAL_CUBE (r * REAL_RECIP (damp));
      if (damp > -50) {
         real expdamp = REAL_EXP (damp);
         scale3 = 1 - expdamp;
         scale5 = 1 - expdamp * (1 - damp);
      }
   }
}

#pragma acc routine seq
CUDA_DEVICE_FUNCTION
inline void damp_thole3 (real r, real pdi, real pti, real pdk, real ptk,
                         real& RESTRICT scale3, real& RESTRICT scale5,
                         real& RESTRICT scale7)
{
   scale3 = 1;
   scale5 = 1;
   scale7 = 1;
   real damp = pdi * pdk;
   if (damp != 0) {
      real pgamma = REAL_MIN (pti, ptk);
      damp = -pgamma * REAL_CUBE (r * REAL_RECIP (damp));
      if (damp > -50) {
         real expdamp = REAL_EXP (damp);
         scale3 = 1 - expdamp;
         scale5 = 1 - expdamp * (1 - damp);
         scale7 = 1 - expdamp * (1 - damp + (real)0.6 * REAL_SQ (damp));
      }
   }
}

#pragma acc routine seq
CUDA_DEVICE_FUNCTION
inline void damp_thole3g (real r, real rr2, real xr, real yr, real zr, real pdi,
                          real pti, real pdk, real ptk, real& RESTRICT scale31,
                          real& RESTRICT scale51, real& RESTRICT scale71,
                          real& RESTRICT rc31, real& RESTRICT rc32,
                          real& RESTRICT rc33, real& RESTRICT rc51,
                          real& RESTRICT rc52, real& RESTRICT rc53,
                          real& RESTRICT rc71, real& RESTRICT rc72,
                          real& RESTRICT rc73)
{
   scale31 = 0;
   scale51 = 0;
   scale71 = 0;
   rc31 = 0;
   rc32 = 0;
   rc33 = 0;
   rc51 = 0;
   rc52 = 0;
   rc53 = 0;
   rc71 = 0;
   rc72 = 0;
   rc73 = 0;
   real damp = pdi * pdk;
   if (damp != 0) {
      real pgamma = REAL_MIN (pti, ptk);
      damp = pgamma * REAL_CUBE (r * REAL_RECIP (damp));
      if (damp < 50) {
         scale31 = REAL_EXP (-damp);
         scale51 = scale31 * (1 + damp);
         scale71 = scale31 * (1 + damp + (real)0.6 * REAL_SQ (damp));
         real temp3 = 3 * damp * scale31 * rr2;
         real temp5 = damp;
         real temp7 = (real)-0.2 + (real)0.6 * damp;
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

#pragma acc routine seq
CUDA_DEVICE_FUNCTION
inline void damp_ewald (real* RESTRICT bn, int order, real r, real rinv,
                        real rr2, real aewald)
{
   real ralpha = aewald * r;
   bn[0] = REAL_ERFC (ralpha) * rinv;
   real alsq2 = 2 * REAL_SQ (aewald);
   real alsq2n = (aewald > 0 ? REAL_RECIP (sqrtpi * aewald) : 0);
   real exp2a = REAL_EXP (-REAL_SQ (ralpha));
   for (int j = 1; j < order; ++j) {
      alsq2n *= alsq2;
      bn[j] = ((j + j - 1) * bn[j - 1] + alsq2n * exp2a) * rr2;
   }
}
TINKER_NAMESPACE_END
