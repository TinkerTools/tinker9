#pragma once
#include "mathfunc.h"
#include "seq_def.h"


namespace tinker {
#pragma acc routine seq
template <int order>
SEQ_CUDA
inline void damp_rep(real* restrict dmpik, real r, real rinv, real r2, real rr3,
                     real rr5, real rr7, real rr9, real rr11, real dmpi,
                     real dmpk)
{
   real r3 = r2 * r;
   real r4 = r2 * r2;
   real r5 = r3 * r2;

   // This makes sure that dmpi > dmpk, which solves numerical issues
   // with atomi,atomk vs. atomk,atomi.
   real dmpsmall = REAL_MIN(dmpi, dmpk);
   real dmpbig = REAL_MAX(dmpi, dmpk);
   dmpi = dmpsmall;
   dmpk = dmpbig;

   real dmpi2 = 0.5f * dmpi;
   real dampi = dmpi2 * r;
   real expi = REAL_EXP(-dampi);
   real dmpi22 = dmpi2 * dmpi2;
   real dmpi23 = dmpi22 * dmpi2;
   real dmpi24 = dmpi22 * dmpi22;
   real dmpi25 = dmpi23 * dmpi22;

   real diff = REAL_ABS(dmpi - dmpk);
   real eps = 0.005f;

   // divisions
   const real div3 = 1 / ((real)3);
   const real div9 = 1 / ((real)9);
   const real div945 = 1 / ((real)945);

   real pre, s, ds, d2s, d3s, d4s, d5s;
   if (diff < eps) {
      real r6 = r3 * r3;
      real r7 = r4 * r3;

      real dmpi26 = dmpi23 * dmpi23;

      pre = 128;
      s = (r + dmpi2 * r2 + dmpi22 * r3 * div3) * expi;
      ds = (dmpi22 * r3 + dmpi23 * r4) * expi * div3;
      d2s = dmpi24 * expi * r5 * div9;
      d3s = dmpi25 * expi * r6 / ((real)45);
      d4s = (dmpi25 * r6 + dmpi26 * r7) * expi / ((real)315);

      if CONSTEXPR (order > 9) {
         real r8 = r4 * r4;
         real dmpi27 = dmpi24 * dmpi23;
         d5s = (dmpi25 * r6 + dmpi26 * r7 + dmpi27 * r8 * div3) * expi * div945;
      }
   } else {
      // treat the case where alpha damping exponents are unequal

      // divisions
      real div5 = 1 / ((real)5);
      real div7 = 1 / ((real)7);
      real div15 = 1 / ((real)15);
      real div21 = 1 / ((real)21);
      real div63 = 1 / ((real)63);
      real div105 = 1 / ((real)105);
      real div189 = 1 / ((real)189);

      real dmpk2 = 0.5f * dmpk;
      real dampk = dmpk2 * r;
      real expk = REAL_EXP(-dampk);
      real dmpk22 = dmpk2 * dmpk2;
      real dmpk23 = dmpk22 * dmpk2;
      real dmpk24 = dmpk22 * dmpk22;
      real dmpk25 = dmpk23 * dmpk22;

      real term = 0.25f * (dmpi + dmpk) * (dmpi - dmpk);
      real term1 = (dmpi + dmpk) * (dmpi - dmpk);
      real tmp = 4 * (dmpi * dmpk) / term1;
      pre = (8192 * dmpi23 * dmpk23) / (term * term * term * term);

      real coef1 = 16  / term1;

      s = (dampi * expk) + (dampk * expi) + tmp * (expi - expk);
      ds =
         (term * dmpk2 * r2 - 4 * (dmpk22 * r + dmpk2)) * dmpi2 * expk / term +
         (term * dmpi2 * r2 + 4 * (dmpi22 * r + dmpi2)) * dmpk2 * expi / term;

      d2s = ((dmpk2 * r2 + dmpk22 * r3) * div3 -
             coef1 * (dmpk23 * r2 * div3 + dmpk22 * r + dmpk2)) *
            dmpi2 * expk +
         ((dmpi2 * r2 + dmpi22 * r3) * div3 +
          coef1 * (dmpi23 * r2 * div3 + dmpi22 * r + dmpi2)) *
            dmpk2 * expi;

      d3s = ((dmpk23 * r4 * div3 + dmpk22 * r3 + dmpk2 * r2) * div5 -
             coef1 *
                (dmpk24 * r3 * div15 + (2 * div5) * dmpk23 * r2 + dmpk22 * r +
                 dmpk2)) *
            dmpi2 * expk +
         ((dmpi23 * r4 * div3 + dmpi22 * r3 + dmpi2 * r2) * div5 +
          coef1 *
             (dmpi24 * r3 * div15 + (2 * div5) * dmpi23 * r2 + dmpi22 * r +
              dmpi2)) *
            dmpk2 * expi;

      d4s = ((dmpk24 * r5 * div15 + 2 * div5 * dmpk23 * r4 + dmpk22 * r3 +
              dmpk2 * r2) *
                div7 -
             coef1 *
                (dmpk25 * r4 * div105 + 2 * div21 * dmpk24 * r3 +
                 3 * div7 * dmpk23 * r2 + dmpk22 * r + dmpk2)) *
            dmpi2 * expk +
         ((dmpi24 * r5 * div15 + 2 * div5 * dmpi23 * r4 + dmpi22 * r3 +
           dmpi2 * r2) *
             div7 +
          coef1 *
             (dmpi25 * r4 * div105 + 2 * div21 * dmpi24 * r3 +
              3 * div7 * dmpi23 * r2 + dmpi22 * r + dmpi2)) *
            dmpk2 * expi;

      if CONSTEXPR (order > 9) {
         real r6 = r3 * r3;
         real dmpi26 = dmpi23 * dmpi23;
         real dmpk26 = dmpk23 * dmpk23;
         d5s = (dmpk25 * r6 * div945 + 2 * div189 * dmpk24 * r5 +
                dmpk23 * r4 * div21 + dmpk22 * r3 * div9 + dmpk2 * r2 * div9 -
                coef1 *
                   (dmpk26 * r5 * div945 + dmpk25 * r4 * div63 +
                    dmpk24 * r3 * div9 + 4 * div9 * dmpk23 * r2 + dmpk22 * r +
                    dmpk2)) *
               dmpi2 * expk +
            (dmpi25 * r6 * div945 + 2 * div189 * dmpi24 * r5 +
             dmpi23 * r4 * div21 + dmpi22 * r3 * div9 + dmpi2 * r2 * div9 +
             coef1 *
                (dmpi26 * r5 * div945 + dmpi25 * r4 * div63 +
                 dmpi24 * r3 * div9 + 4 * div9 * dmpi23 * r2 + dmpi22 * r +
                 dmpi2)) *
               dmpk2 * expi;
      }
   }

   // convert partial derivatives into full derivatives
   s *= rinv;
   ds *= rr3;
   d2s *= rr5;
   d3s *= rr7;
   d4s *= rr9;

   dmpik[0] = 0.5f * pre * s * s;
   dmpik[1] = pre * s * ds;
   dmpik[2] = pre * (s * d2s + ds * ds);
   dmpik[3] = pre * (s * d3s + 3 * ds * d2s);
   dmpik[4] = pre * (s * d4s + 4 * ds * d3s + 3 * d2s * d2s);

   if CONSTEXPR (order > 9) {
      d5s *= rr11;
      dmpik[5] = pre * (s * d5s + 5 * ds * d4s + 10 * d2s * d3s);
   }
}
}
