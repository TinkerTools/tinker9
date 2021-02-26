#pragma once
#include "mathfunc.h"
#include "seq_def.h"


namespace tinker {
#pragma acc routine seq
template <int ORDER>
SEQ_CUDA
inline void damp_pole_v1(real* restrict dmpik, real* restrict dmpi,
                         real* restrict dmpk, real r, real alphai, real alphak)
{
#if TINKER_REAL_SIZE == 8
   real eps = 0.001f;
#elif TINKER_REAL_SIZE == 4
   real eps = 0.005f;
#endif

   real diff = REAL_ABS(alphai - alphak);

   if (diff < eps)
      alphai = 0.5f * (alphai + alphak);

   real dampi = alphai * r;
   real dampk = alphak * r;
   real expi = REAL_EXP(-dampi);
   real expk = REAL_EXP(-dampk);

   real dampi2 = dampi * dampi;
   real dampi3 = dampi * dampi2;
   real dampi4 = dampi2 * dampi2;
   real dampi5 = dampi2 * dampi3;

   // divisions
   const real div3 = 1 / ((real)3);
   const real div6 = 1 / ((real)6);
   const real div15 = 1 / ((real)15);
   const real div21 = 1 / ((real)21);
   const real div30 = 1 / ((real)30);
   const real div105 = 1 / ((real)105);

   // GORDON1
   // core-valence
   dmpi[0] = 1 - (1 + 0.5f * dampi) * expi;
   dmpi[1] = 1 - (1 + dampi + 0.5f * dampi2) * expi;
   dmpi[2] = 1 - (1 + dampi + 0.5f * dampi2 + dampi3 * div6) * expi;
   dmpi[3] =
      1 - (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + dampi4 * div30) * expi;
   dmpi[4] = 1 -
      (1 + dampi + 0.5f * dampi2 + dampi3 * div6 +
       (4 * dampi4 + 0.5f * dampi5) * div105) *
         expi;

   if (diff < eps) {
      dmpk[0] = dmpi[0];
      dmpk[1] = dmpi[1];
      dmpk[2] = dmpi[2];
      dmpk[3] = dmpi[3];
      dmpk[4] = dmpi[4];

      // valence-valence
      real dampi6 = dampi3 * dampi3;
      real dampi7 = dampi3 * dampi4;
      const real div5 = 1 / ((real)5);
      const real div16 = 1 / ((real)16);
      const real div24 = 1 / ((real)24);
      const real div42 = 1 / ((real)42);
      const real div48 = 1 / ((real)48);
      const real div120 = 1 / ((real)120);
      const real div144 = 1 / ((real)144);

      dmpik[0] =
         1 - (1 + (11 * dampi + 3 * dampi2) * div16 + dampi3 * div48) * expi;
      dmpik[1] =
         1 - (1 + dampi + 0.5f * dampi2 + (7 * dampi3 + dampi4) * div48) * expi;
      dmpik[2] = 1 -
         (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + dampi4 * div24 +
          dampi5 * div144) *
            expi;
      dmpik[3] = 1 -
         (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + dampi4 * div24 +
          (dampi5 + dampi6 * div6) * div120) *
            expi;
      dmpik[4] = 1 -
         ((1 + dampi + 0.5f * dampi2 + dampi3 * div6) +
          ((dampi4 + dampi5 * div5) +
           0.1f * (dampi6 * div3 + dampi7 * div21)) *
             div24) *
            expi;
      if CONSTEXPR (ORDER > 9) {
         real dampi8 = dampi4 * dampi4;
         const real div378 = 1 / ((real)378);
         dmpik[5] = 1 -
            (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + dampi4 * div24 +
             (dampi5 + dampi6 * div6 + dampi7 * div42 + dampi8 * div378) *
                div120) *
               expi;
      }
   } else {
      real dampk2 = dampk * dampk;
      real dampk3 = dampk * dampk2;
      real dampk4 = dampk2 * dampk2;
      real dampk5 = dampk2 * dampk3;

      const real div5 = 1 / ((real)5);
      const real div7 = 1 / ((real)7);

      dmpk[0] = 1 - (1 + 0.5f * dampk) * expk;
      dmpk[1] = 1 - (1 + dampk + 0.5f * dampk2) * expk;
      dmpk[2] = 1 - (1 + dampk + 0.5f * dampk2 + dampk3 * div6) * expk;
      dmpk[3] = 1 -
         (1 + dampk + 0.5f * dampk2 + dampk3 * div6 + dampk4 * div30) * expk;
      dmpk[4] = 1 -
         (1 + dampk + 0.5f * dampk2 + dampk3 * div6 +
          (4 * dampk4 + 0.5f * dampk5) * div105) *
            expk;

      // valence-valence
      real alphai2 = alphai * alphai;
      real alphak2 = alphak * alphak;
      real alphaik = ((alphak + alphai) * (alphak - alphai));
      real termi = alphak2 / alphaik;
      real termk = -alphai2 / alphaik;
      real termi2 = termi * termi;
      real termk2 = termk * termk;

      dmpik[0] = 1 - termi2 * (1 + 2 * termk + 0.5f * dampi) * expi -
         termk2 * (1 + 2 * termi + 0.5f * dampk) * expk;
      dmpik[1] = 1 - termi2 * (1 + dampi + 0.5f * dampi2) * expi -
         termk2 * (1 + dampk + 0.5f * dampk2) * expk -
         2 * termi2 * termk * (1 + dampi) * expi -
         2 * termk2 * termi * (1 + dampk) * expk;
      dmpik[2] = 1 -
         termi2 * (1 + dampi + 0.5f * dampi2 + dampi3 * div6) * expi -
         termk2 * (1 + dampk + 0.5f * dampk2 + dampk3 * div6) * expk -
         2 * termi2 * termk * (1 + dampi + dampi2 * div3) * expi -
         2 * termk2 * termi * (1 + dampk + dampk2 * div3) * expk;
      dmpik[3] = 1 -
         termi2 * (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + dampi4 * div30) *
            expi -
         termk2 * (1 + dampk + 0.5f * dampk2 + dampk3 * div6 + dampk4 * div30) *
            expk -
         2 * termi2 * termk * (1 + dampi + 2 * dampi2 * div5 + dampi3 * div15) *
            expi -
         2 * termk2 * termi * (1 + dampk + 2 * dampk2 * div5 + dampk3 * div15) *
            expk;
      dmpik[4] = 1 -
         termi2 *
            (1 + dampi + 0.5f * dampi2 + dampi3 * div6 +
             (4 * dampi4 + 0.5f * dampi5) * div105) *
            expi -
         termk2 *
            (1 + dampk + 0.5f * dampk2 + dampk3 * div6 +
             (4 * dampk4 + 0.5f * dampk5) * div105) *
            expk -
         2 * termi2 * termk *
            (1 + dampi + 3 * dampi2 * div7 + 2 * dampi3 * div21 +
             dampi4 * div105) *
            expi -
         2 * termk2 * termi *
            (1 + dampk + 3 * dampk2 * div7 + 2 * dampk3 * div21 +
             dampk4 * div105) *
            expk;

      if CONSTEXPR (ORDER > 9) {
         real dampi6 = dampi3 * dampi3;
         real dampk6 = dampk3 * dampk3;
         const real div945 = 1 / ((real)945);
         const real div9 = 1 / ((real)9);
         const real div63 = 1 / ((real)63);
         const real div126 = 1 / ((real)126);
         const real div315 = 1 / ((real)315);
         const real div1890 = 1 / ((real)1890);

         dmpik[5] = 1 -
            (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + 5 * dampi4 * div126 +
             2 * dampi5 * div315 + dampi6 * div1890) *
               termi2 * expi -
            (1 + dampk + 0.5f * dampk2 + dampk3 * div6 + 5 * dampk4 * div126 +
             2 * dampk5 * div315 + dampk6 * div1890) *
               termk2 * expk -
            (1 + dampi + (4 * dampi2 + dampi3) * div9 + dampi4 * div63 +
             dampi5 * div945) *
               2 * termi2 * termk * expi -
            (1 + dampk + (4 * dampk2 + dampk3) * div9 + dampk4 * div63 +
             dampk5 * div945) *
               2 * termk2 * termi * expk;
      }
   }
}


#pragma acc routine seq
template <int ORDER>
SEQ_CUDA
inline void damp_pole_v2(real* restrict dmpik, real* restrict dmpi,
                         real* restrict dmpk, real r, real alphai, real alphak)
{
#if TINKER_REAL_SIZE == 8
   real eps = 0.001f;
#elif TINKER_REAL_SIZE == 4
   real eps = 0.05f;
#endif

   real diff = REAL_ABS(alphai - alphak);

   if (diff < eps)
      alphai = 0.5f * (alphai + alphak);

   real dampi = alphai * r;
   real dampk = alphak * r;
   real expi = REAL_EXP(-dampi);
   real expk = REAL_EXP(-dampk);

   real dampi2 = dampi * dampi;
   real dampi3 = dampi * dampi2;
   real dampi4 = dampi2 * dampi2;
   real dampi5 = dampi2 * dampi3;

   // divisions
   const real div3 = 1 / ((real)3);
   const real div6 = 1 / ((real)6);
   const real div15 = 1 / ((real)15);
   const real div21 = 1 / ((real)21);
   const real div30 = 1 / ((real)30);
   const real div105 = 1 / ((real)105);

   // GORDON1
   // core-valence
   dmpi[0] = (1 + 0.5f * dampi) * expi;
   dmpi[1] = (1 + dampi + 0.5f * dampi2) * expi;
   dmpi[2] = (1 + dampi + 0.5f * dampi2 + dampi3 * div6) * expi;
   dmpi[3] =
      (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + dampi4 * div30) * expi;
   dmpi[4] = (1 + dampi + 0.5f * dampi2 + dampi3 * div6 +
              (4 * dampi4 + 0.5f * dampi5) * div105) *
      expi;

   if (diff < eps) {
      dmpk[0] = dmpi[0];
      dmpk[1] = dmpi[1];
      dmpk[2] = dmpi[2];
      dmpk[3] = dmpi[3];
      dmpk[4] = dmpi[4];

      // valence-valence
      real dampi6 = dampi3 * dampi3;
      real dampi7 = dampi3 * dampi4;
      const real div5 = 1 / ((real)5);
      const real div16 = 1 / ((real)16);
      const real div24 = 1 / ((real)24);
      const real div42 = 1 / ((real)42);
      const real div48 = 1 / ((real)48);
      const real div120 = 1 / ((real)120);
      const real div144 = 1 / ((real)144);

      dmpik[0] =
         (1 + (11 * dampi + 3 * dampi2) * div16 + dampi3 * div48) * expi;
      dmpik[1] =
         (1 + dampi + 0.5f * dampi2 + (7 * dampi3 + dampi4) * div48) * expi;
      dmpik[2] = (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + dampi4 * div24 +
                  dampi5 * div144) *
         expi;
      dmpik[3] = (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + dampi4 * div24 +
                  (dampi5 + dampi6 * div6) * div120) *
         expi;
      dmpik[4] =
         ((1 + dampi + 0.5f * dampi2 + dampi3 * div6) +
          ((dampi4 + dampi5 * div5) + 0.1f * (dampi6 * div3 + dampi7 * div21)) *
             div24) *
         expi;
      if CONSTEXPR (ORDER > 9) {
         real dampi8 = dampi4 * dampi4;
         const real div378 = 1 / ((real)378);
         dmpik[5] =
            (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + dampi4 * div24 +
             (dampi5 + dampi6 * div6 + dampi7 * div42 + dampi8 * div378) *
                div120) *
            expi;
      }
   } else {
      real dampk2 = dampk * dampk;
      real dampk3 = dampk * dampk2;
      real dampk4 = dampk2 * dampk2;
      real dampk5 = dampk2 * dampk3;

      const real div5 = 1 / ((real)5);
      const real div7 = 1 / ((real)7);

      dmpk[0] = (1 + 0.5f * dampk) * expk;
      dmpk[1] = (1 + dampk + 0.5f * dampk2) * expk;
      dmpk[2] = (1 + dampk + 0.5f * dampk2 + dampk3 * div6) * expk;
      dmpk[3] =
         (1 + dampk + 0.5f * dampk2 + dampk3 * div6 + dampk4 * div30) * expk;
      dmpk[4] = (1 + dampk + 0.5f * dampk2 + dampk3 * div6 +
                 (4 * dampk4 + 0.5f * dampk5) * div105) *
         expk;

      // valence-valence
      real alphai2 = alphai * alphai;
      real alphak2 = alphak * alphak;
      real alphaik = ((alphak + alphai) * (alphak - alphai));
      real termi = alphak2 / alphaik;
      real termk = -alphai2 / alphaik;
      real termi2 = termi * termi;
      real termk2 = termk * termk;

      dmpik[0] = termi2 * (1 + 2 * termk + 0.5f * dampi) * expi +
         termk2 * (1 + 2 * termi + 0.5f * dampk) * expk;
      dmpik[1] = termi2 * (1 + dampi + 0.5f * dampi2) * expi +
         termk2 * (1 + dampk + 0.5f * dampk2) * expk +
         2 * termi2 * termk * (1 + dampi) * expi +
         2 * termk2 * termi * (1 + dampk) * expk;
      dmpik[2] = termi2 * (1 + dampi + 0.5f * dampi2 + dampi3 * div6) * expi +
         termk2 * (1 + dampk + 0.5f * dampk2 + dampk3 * div6) * expk +
         2 * termi2 * termk * (1 + dampi + dampi2 * div3) * expi +
         2 * termk2 * termi * (1 + dampk + dampk2 * div3) * expk;
      dmpik[3] = termi2 *
            (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + dampi4 * div30) *
            expi +
         termk2 * (1 + dampk + 0.5f * dampk2 + dampk3 * div6 + dampk4 * div30) *
            expk +
         2 * termi2 * termk * (1 + dampi + 2 * dampi2 * div5 + dampi3 * div15) *
            expi +
         2 * termk2 * termi * (1 + dampk + 2 * dampk2 * div5 + dampk3 * div15) *
            expk;
      dmpik[4] = termi2 *
            (1 + dampi + 0.5f * dampi2 + dampi3 * div6 +
             (4 * dampi4 + 0.5f * dampi5) * div105) *
            expi +
         termk2 *
            (1 + dampk + 0.5f * dampk2 + dampk3 * div6 +
             (4 * dampk4 + 0.5f * dampk5) * div105) *
            expk +
         2 * termi2 * termk *
            (1 + dampi + 3 * dampi2 * div7 + 2 * dampi3 * div21 +
             dampi4 * div105) *
            expi +
         2 * termk2 * termi *
            (1 + dampk + 3 * dampk2 * div7 + 2 * dampk3 * div21 +
             dampk4 * div105) *
            expk;

      if CONSTEXPR (ORDER > 9) {
         real dampi6 = dampi3 * dampi3;
         real dampk6 = dampk3 * dampk3;
         const real div945 = 1 / ((real)945);
         const real div9 = 1 / ((real)9);
         const real div63 = 1 / ((real)63);
         const real div126 = 1 / ((real)126);
         const real div315 = 1 / ((real)315);
         const real div1890 = 1 / ((real)1890);

         dmpik[5] =
            (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + 5 * dampi4 * div126 +
             2 * dampi5 * div315 + dampi6 * div1890) *
               termi2 * expi +
            (1 + dampk + 0.5f * dampk2 + dampk3 * div6 + 5 * dampk4 * div126 +
             2 * dampk5 * div315 + dampk6 * div1890) *
               termk2 * expk +
            (1 + dampi + (4 * dampi2 + dampi3) * div9 + dampi4 * div63 +
             dampi5 * div945) *
               2 * termi2 * termk * expi +
            (1 + dampk + (4 * dampk2 + dampk3) * div9 + dampk4 * div63 +
             dampk5 * div945) *
               2 * termk2 * termi * expk;
      }
   }
}


#pragma acc routine seq
SEQ_CUDA
inline void damp_dir(real* restrict dmpi, real* restrict dmpk, real r,
                     real alphai, real alphak)
{
   real eps = 0.001f;
   real diff = REAL_ABS(alphai - alphak);
   real dampi = alphai * r;
   real dampk = alphak * r;
   real expi = REAL_EXP(-dampi);
   real expk = REAL_EXP(-dampk);

   real dampi2 = dampi * dampi;
   real dampi3 = dampi * dampi2;

   // GORDON1
   real dampi4 = dampi2 * dampi2;

   // divisions
   const real div6 = 1 / ((real)6);
   const real div30 = 1 / ((real)30);

   // core-valence
   dmpi[1] = 1 - (1 + dampi + 0.5f * dampi2) * expi;
   dmpi[2] = 1 - (1 + dampi + 0.5f * dampi2 + dampi3 * div6) * expi;
   dmpi[3] =
      1 - (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + dampi4 * div30) * expi;

   if (diff < eps) {
      dmpk[1] = dmpi[1];
      dmpk[2] = dmpi[2];
      dmpk[3] = dmpi[3];
   } else {
      real dampk2 = dampk * dampk;
      real dampk3 = dampk * dampk2;
      real dampk4 = dampk2 * dampk2;
      dmpk[1] = 1 - (1 + dampk + 0.5f * dampk2) * expk;
      dmpk[2] = 1 - (1 + dampk + 0.5f * dampk2 + dampk3 * div6) * expk;
      dmpk[3] = 1 -
         (1 + dampk + 0.5f * dampk2 + dampk3 * div6 + dampk4 * div30) * expk;
   }
}


#pragma acc routine seq
SEQ_CUDA
inline void damp_mut(real* restrict dmpik, real r, real alphai, real alphak)
{
#if TINKER_REAL_SIZE == 8
   real eps = 0.001f;
#elif TINKER_REAL_SIZE == 4
   real eps = 0.05f;
#endif

   real diff = REAL_ABS(alphai - alphak);

   if (diff < eps)
      alphai = 0.5f * (alphai + alphak);

   real dampi = alphai * r;
   real dampk = alphak * r;
   real expi = REAL_EXP(-dampi);
   real expk = REAL_EXP(-dampk);

   real dampi2 = dampi * dampi;

   // divisions
   const real div6 = 1 / ((real)6);

   // GORDON1
   real dampi3 = dampi * dampi2;

   if (diff < eps) {
      real dampi4 = dampi2 * dampi2;
      real dampi5 = dampi2 * dampi3;
      const real div24 = 1 / ((real)24);
      const real div48 = 1 / ((real)48);
      const real div144 = 1 / ((real)144);

      dmpik[1] =
         1 - (1 + dampi + 0.5f * dampi2 + (7 * dampi3 + dampi4) * div48) * expi;
      dmpik[2] = 1 -
         (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + dampi4 * div24 +
          dampi5 * div144) *
            expi;
   } else {
      real dampk2 = dampk * dampk;
      real dampk3 = dampk * dampk2;
      real alphai2 = alphai * alphai;
      real alphak2 = alphak * alphak;
      real termi = alphak2 / (alphak2 - alphai2);
      real termk = alphai2 / (alphai2 - alphak2);
      real termi2 = termi * termi;
      real termk2 = termk * termk;

      const real div3 = 1 / ((real)3);

      dmpik[1] = 1 - termi2 * (1 + dampi + 0.5f * dampi2) * expi -
         termk2 * (1 + dampk + 0.5f * dampk2) * expk -
         2 * termi2 * termk * (1 + dampi) * expi -
         2 * termk2 * termi * (1 + dampk) * expk;
      dmpik[2] = 1 -
         termi2 * (1 + dampi + 0.5f * dampi2 + dampi3 * div6) * expi -
         termk2 * (1 + dampk + 0.5f * dampk2 + dampk3 * div6) * expk -
         2 * termi2 * termk * (1 + dampi + dampi2 * div3) * expi -
         2 * termk2 * termi * (1 + dampk + dampk2 * div3) * expk;
   }
}


#pragma acc routine seq
SEQ_CUDA
inline void damp_pot(real* restrict dmpk, real r, real alphak)
{
   real dampk = alphak * r;
   real expk = REAL_EXP(-dampk);
   real dampk2 = dampk * dampk;

   // GORDON1
   // divisions
   const real div6 = 1 / ((real)6);

   real dampk3 = dampk * dampk2;
   dmpk[0] = 1 - (1 + 0.5f * dampk) * expk;
   dmpk[1] = 1 - (1 + dampk + 0.5f * dampk2) * expk;
   dmpk[2] = 1 - (1 + dampk + 0.5f * dampk2 + dampk3 * div6) * expk;
}


#pragma acc routine seq
template <int order>
SEQ_CUDA
inline void damp_rep(real* restrict dmpik, real r, real rinv, real r2, real rr3,
                     real rr5, real rr7, real rr9, real rr11, real dmpi,
                     real dmpk)
{
#if TINKER_REAL_SIZE == 8
   real eps = 0.001f;
#elif TINKER_REAL_SIZE == 4
   real eps = 0.05f;
#endif

   // This makes sure that dmpi > dmpk, which solves numerical issues
   // with atomi,atomk vs. atomk,atomi.
   real dmpsmall = REAL_MIN(dmpi, dmpk);
   real dmpbig = REAL_MAX(dmpi, dmpk);
   dmpi = dmpsmall;
   dmpk = dmpbig;

   real diff = REAL_ABS(dmpi - dmpk);
   if (diff < eps)
      dmpi = 0.5f * (dmpi + dmpk);

   real r3 = r2 * r;
   real r4 = r2 * r2;
   real r5 = r3 * r2;
   real dmpi2 = 0.5f * dmpi;
   real dampi = dmpi2 * r;
   real expi = REAL_EXP(-dampi);
   real dmpi22 = dmpi2 * dmpi2;
   real dmpi23 = dmpi22 * dmpi2;
   real dmpi24 = dmpi22 * dmpi22;
   real dmpi25 = dmpi23 * dmpi22;


   // divisions
   const real div3 = 1 / ((real)3);
   const real div9 = 1 / ((real)9);
   const real div945 = 1 / ((real)945);

   real pre, s, ds, d2s, d3s, d4s, d5s;
   if (diff < eps) {
      real r6 = r3 * r3;
      real r7 = r4 * r3;

      real dmpi26 = dmpi23 * dmpi23;

      pre = 2;
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
      pre = (128 * dmpi23 * dmpk23) / (term * term * term * term);

      real coef1 = 16 / term1;

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
