#pragma once
#include "mathfunc.h"
#include "seq_def.h"


namespace tinker {
#pragma acc routine seq
template <int ORDER>
SEQ_CUDA
inline void damp_pole(real* restrict dmpik, real* restrict dmpi,
                      real* restrict dmpk, real r, real alphai, real alphak)
{
   real eps = 0.001f;
   real diff = REAL_ABS(alphai - alphak);
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
   const real div30 = 1 / ((real)30);
   const real div105 = 1 / ((real)105);
   const real div210 = 1 / ((real)210);

   // if GORDON1
   // if CONSTEXPR (eq<PENTYP, GORDON1>()) {
   // core-valence
   dmpi[0] = 1 - (1 + 0.5f * dampi) * expi;
   dmpi[1] = 1 - (1 + dampi + 0.5f * dampi2) * expi;
   dmpi[2] = 1 - (1 + dampi + 0.5f * dampi2 + dampi3 * div6) * expi;
   dmpi[3] =
      1 - (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + dampi4 * div30) * expi;
   dmpi[4] = 1 -
      (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + 4 * dampi4 * div105 +
       dampi5 * div210) *
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
      const real div16 = 1 / ((real)16);
      const real div24 = 1 / ((real)24);
      const real div48 = 1 / ((real)48);
      const real div120 = 1 / ((real)120);
      const real div144 = 1 / ((real)144);
      const real div720 = 1 / ((real)720);
      const real div5040 = 1 / ((real)5040);


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
          dampi5 * div120 + dampi6 * div720) *
            expi;
      dmpik[4] = 1 -
         (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + dampi4 * div24 +
          dampi5 * div120 + dampi6 * div720 + dampi7 * div5040) *
            expi;

      if CONSTEXPR (ORDER > 9) {
         real dampi8 = dampi4 * dampi4;
         const real div45360 = 1 / ((real)45360);
         dmpik[5] = 1 -
            (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + dampi4 * div24 +
             dampi5 * div120 + dampi6 * div720 + dampi7 * div5040 +
             div45360 * dampi8) *
               expi;
      }

   } else {
      real dampk2 = dampk * dampk;
      real dampk3 = dampk * dampk2;
      real dampk4 = dampk2 * dampk2;
      real dampk5 = dampk2 * dampk3;

      const real div5 = 1 / ((real)5);
      const real div7 = 1 / ((real)7);
      const real div21 = 1 / ((real)21);
      const real div210 = 1 / ((real)210);

      dmpk[0] = 1 - (1 + 0.5f * dampk) * expk;
      dmpk[1] = 1 - (1 + dampk + 0.5f * dampk2) * expk;
      dmpk[2] = 1 - (1 + dampk + 0.5f * dampk2 + dampk3 * div6) * expk;
      dmpk[3] = 1 -
         (1 + dampk + 0.5f * dampk2 + dampk3 * div6 + dampk4 * div30) * expk;
      dmpk[4] = 1 -
         (1 + dampk + 0.5f * dampk2 + +dampk3 * div6 + 4 * dampk4 * div105 +
          dampk5 * div210) *
            expk;


      // valence-valence
      real alphai2 = alphai * alphai;
      real alphak2 = alphak * alphak;
      real termi = alphak2 / (alphak2 - alphai2);
      real termk = alphai2 / (alphai2 - alphak2);
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
            (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + 4 * dampi4 * div105 +
             dampi5 * div210) *
            expi -
         termk2 *
            (1 + dampk + 0.5f * dampk2 + dampk3 * div6 + 4 * dampk4 * div105 +
             dampk5 * div210) *
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
             2 * dampi5 * div315 + dampi6 * div1890) * termi2 * expi -
            (1 + dampk + 0.5f * dampk2 + dampk3 * div6 + 5 * dampk4 * div126 +
             2 * dampk5 * div315 + dampk6 * div1890) * termk2 * expk -
            (1 + dampi + (4 * dampi2 + dampi3) * div9 + dampi4 * div63 +
               dampi5 * div945) * 2 * termi2 * termk * expi -
            (1 + dampk + (4 * dampk2 + dampk3) * div9 + dampk4 * div63 +
               dampk5 * div945) * 2 * termk2 * termi * expk;
      }
   }
   //} // endif GORDON1


   // if CONSTEXPR (eq<PENTYP, GORDON2>()) {
   //    // core-valence
   //    dmpi[0] = 1 - expi;
   //    dmpi[1] = 1 - (1 + dampi) * expi;
   //    dmpi[2] = 1 - (1 + dampi + dampi2 * div3) * expi;
   //    dmpi[3] = 1 - (1 + dampi + 0.4f * dampi2 + dampi3 * div15) * expi;

   //    if (diff < eps) {
   //       dmpk[0] = dmpi[0];
   //       dmpk[1] = dmpi[1];
   //       dmpk[2] = dmpi[2];
   //       dmpk[3] = dmpi[3];

   //       // valence-valence
   //       real dampi6 = dampi3 * dampi3;
   //       const real div210 = 1 / ((real)210);

   //       dmpik[0] = 1 - (1 + 0.5f * dampi) * expi;
   //       dmpik[1] = 1 - (1 + dampi + 0.5f * dampi2) * expi;
   //       dmpik[2] = 1 - (1 + dampi + 0.5f * dampi2 + dampi3 * div6) * expi;
   //       dmpik[3] = 1 -
   //          (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + dampi4 * div30) *
   //          expi;
   //       dmpik[4] = 1 -
   //          (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + 4 * dampi4 * div105
   //          +
   //           dampi5 * div210) *
   //             expi;

   //       if CONSTEXPR (ORDER > 9) {
   //          const real div126 = 1 / ((real)126);
   //          const real div315 = 1 / ((real)315);
   //          const real div1890 = 1 / ((real)1890);
   //          dmpik[5] = 1 -
   //             (1 + dampi + 0.5f * dampi2 + dampi3 * div6 +
   //              5 * dampi4 * div126 + 2 * dampi5 * div315 + dampi6 * div1890)
   //              *
   //                expi;
   //       }


   //    } else {
   //       real dampk2 = dampk * dampk;
   //       real dampk3 = dampk * dampk2;
   //       real dampk4 = dampk2 * dampk2;
   //       real dampk5 = dampk2 * dampk3;

   //       const real div7 = 1 / ((real)7);
   //       const real div21 = 1 / ((real)21);

   //       // core-valence
   //       dmpk[0] = 1 - expk;
   //       dmpk[1] = 1 - (1 + dampk) * expk;
   //       dmpk[2] = 1 - (1 + dampk + dampk2 * div3) * expk;
   //       dmpk[3] = 1 - (1 + dampk + 0.4f * dampk2 + dampk3 * div15) * expk;

   //       // valence-valence
   //       real alphai2 = alphai * alphai;
   //       real alphak2 = alphak * alphak;
   //       real termi = alphak2 / (alphak2 - alphai2);
   //       real termk = alphai2 / (alphai2 - alphak2);
   //       dmpik[0] = 1 - termi * expi - termk * expk;
   //       dmpik[1] = 1 - termi * (1 + dampi) * expi - termk * (1 + dampk) *
   //       expk; dmpik[2] = 1 - termi * (1 + dampi + dampi2 * div3) * expi -
   //          termk * (1 + dampk + dampk2 * div3) * expk;
   //       dmpik[3] = 1 -
   //          termi * (1 + dampi + 0.4f * dampi2 + dampi3 * div15) * expi -
   //          termk * (1 + dampk + 0.4f * dampk2 + dampk3 * div15) * expk;
   //       dmpik[4] = 1 -
   //          termi *
   //             (1 + dampi + 3 * dampi2 * div7 + 2 * dampi3 * div21 +
   //              dampi4 * div105) *
   //             expi -
   //          termk *
   //             (1 + dampk + 3 * dampk2 * div7 + 2 * dampk3 * div21 +
   //              dampk4 * div105) *
   //             expk;

   //       if CONSTEXPR (ORDER > 9) {
   //          const real div945 = 1 / ((real)945);
   //          const real div9 = 1 / ((real)9);
   //          const real div63 = 1 / ((real)63);

   //          dmpik[5] = 1 -
   //             termi *
   //                (1 + dampi + (4 * dampi2 + dampi3) * div9 + dampi4 * div63
   //                +
   //                 dampi5 * div945) *
   //                expi -
   //             termk *
   //                (1 + dampk + (4 * dampk2 + dampk3) * div9 + dampk4 * div63
   //                +
   //                 dampk5 * div945) *
   //                expk;
   //       }
   //    }

   // } // endif GORDON2
}

#pragma acc routine seq
SEQ_CUDA
inline void damp_polar(real* restrict dmpik, real* restrict dmpi,
                       real* restrict dmpk, real r, real alphai, real alphak)
{
   real eps = 0.001f;
   real diff = REAL_ABS(alphai - alphak);
   real dampi = alphai * r;
   real dampk = alphak * r;
   real expi = REAL_EXP(-dampi);
   real expk = REAL_EXP(-dampk);

   real dampi2 = dampi * dampi;
   real dampi3 = dampi * dampi2;

   // divisions
   const real div3 = 1 / ((real)3);
   const real div6 = 1 / ((real)6);

   // if CONSTEXPR (eq<PENTYP, GORDON1>()) {
   // core-valence
   dmpi[1] = 1 - (1 + dampi + 0.5f * dampi2) * expi;
   dmpi[2] = 1 - (1 + dampi + 0.5f * dampi2 + dampi3 * div6) * expi;

   if (diff < eps) {
      real dampi4 = dampi2 * dampi2;
      real dampi5 = dampi2 * dampi3;
      const real div24 = 1 / ((real)24);
      const real div48 = 1 / ((real)48);
      const real div144 = 1 / ((real)144);

      dmpk[1] = dmpi[1];
      dmpk[2] = dmpi[2];

      // valence-valence
      dmpik[1] =
         1 - (1 + dampi + 0.5f * dampi2 + (7 * dampi3 + dampi4) * div48) * expi;
      dmpik[2] = 1 -
         (1 + dampi + 0.5f * dampi2 + dampi3 * div6 + dampi4 * div24 +
          dampi5 * div144) *
            expi;

   } else {
      real dampk2 = dampk * dampk;
      real dampk3 = dampk * dampk2;
      dmpk[1] = 1 - (1 + dampk + 0.5f * dampk2) * expk;
      dmpk[2] = 1 - (1 + dampk + 0.5f * dampk2 + dampk3 * div6) * expk;

      // valence-valence
      real alphai2 = alphai * alphai;
      real alphak2 = alphak * alphak;
      real termi = alphak2 / (alphak2 - alphai2);
      real termk = alphai2 / (alphai2 - alphak2);
      real termi2 = termi * termi;
      real termk2 = termk * termk;

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

   // } // endif GORDON1


   // if CONSTEXPR (eq<PENTYP, GORDON2>()) {
   //    // core-valence
   //    dmpi[1] = 1 - (1 + dampi) * expi;
   //    dmpi[2] = 1 - (1 + dampi + dampi2 * div3) * expi;

   //    if (diff < eps) {
   //       dmpk[1] = dmpi[1];
   //       dmpk[2] = dmpi[2];

   //       // valence-valence
   //       dmpik[1] = 1 - (1 + dampi + 0.5f * dampi2) * expi;
   //       dmpik[2] = 1 - (1 + dampi + 0.5f * dampi2 + dampi3 * div6) * expi;
   //    } else {
   //       real dampk2 = dampk * dampk;
   //       dmpk[1] = 1 - (1 + dampk) * expk;
   //       dmpk[2] = 1 - (1 + dampk + dampk2 * div3) * expk;

   //       // valence-valence
   //       real alphai2 = alphai * alphai;
   //       real alphak2 = alphak * alphak;
   //       real termi = alphak2 / (alphak2 - alphai2);
   //       real termk = alphai2 / (alphai2 - alphak2);

   //       dmpik[1] = 1 - termi * (1 + dampi) * expi - termk * (1 + dampk) *
   //       expk; dmpik[2] = 1 - termi * (1 + dampi + dampi2 * div3) * expi -
   //          termk * (1 + dampk + dampk2 * div3) * expk;
   //    }

   // } // endif GORDON2
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


   // if CONSTEXPR (eq<PENTYP, GORDON1>()) {
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

   // } // endif GORDON1


   // if CONSTEXPR (eq<PENTYP, GORDON2>()) {
   //    const real div3 = 1 / ((real)3);
   //    const real div15 = 1 / ((real)15);

   //    dmpi[1] = 1 - (1 + dampi) * expi;
   //    dmpi[2] = 1 - (1 + dampi + dampi2 * div3) * expi;
   //    dmpi[3] = 1 - (1 + dampi + 0.4f * dampi2 + dampi3 * div15) * expi;

   //    if (diff < eps) {
   //       dmpk[1] = dmpi[1];
   //       dmpk[2] = dmpi[2];
   //       dmpk[3] = dmpi[3];
   //    } else {
   //       real dampk2 = dampk * dampk;
   //       real dampk3 = dampk * dampk2;
   //       dmpk[1] = 1 - (1 + dampk) * expk;
   //       dmpk[2] = 1 - (1 + dampk + dampk2 * div3) * expk;
   //       dmpk[3] = 1 - (1 + dampk + 0.4f * dampk2 + dampk3 * div15) * expk;
   //    }
   // } // endif GORDON2
}

#pragma acc routine seq
SEQ_CUDA
inline void damp_mut(real* restrict dmpik, real r, real alphai, real alphak)
{
   real eps = 0.001f;
   real diff = REAL_ABS(alphai - alphak);
   real dampi = alphai * r;
   real dampk = alphak * r;
   real expi = REAL_EXP(-dampi);
   real expk = REAL_EXP(-dampk);

   real dampi2 = dampi * dampi;

   // divisions
   const real div6 = 1 / ((real)6);

   // if CONSTEXPR (eq<PENTYP, GORDON1>()) {
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
   // } // endif GORDON1


   // if CONSTEXPR (eq<PENTYP, GORDON2>()) {
   //    if (diff < eps) {
   //       real dampi3 = dampi * dampi2;
   //       dmpik[1] = 1 - (1 + dampi + 0.5f * dampi2) * expi;
   //       dmpik[2] = 1 - (1 + dampi + 0.5f * dampi2 + dampi3 * div6) * expi;
   //    } else {
   //       real dampk2 = dampk * dampk;
   //       real alphai2 = alphai * alphai;
   //       real alphak2 = alphak * alphak;
   //       real termi = alphak2 / (alphak2 - alphai2);
   //       real termk = alphai2 / (alphai2 - alphak2);

   //       const real div3 = 1 / ((real)3);

   //       dmpik[1] = 1 - termi * (1 + dampi) * expi - termk * (1 + dampk) *
   //       expk; dmpik[2] = 1 - termi * (1 + dampi + dampi2 * div3) * expi -
   //          termk * (1 + dampk + dampk2 * div3) * expk;
   //    }
   // } // endif GORDON2
}

#pragma acc routine seq
SEQ_CUDA
inline void damp_pot(real* restrict dmpk, real r, real alphak)
{
   real dampk = alphak * r;
   real expk = REAL_EXP(-dampk);
   real dampk2 = dampk * dampk;


   // if CONSTEXPR (eq<PENTYP, GORDON1>()) {
   // divisions
   const real div6 = 1 / ((real)6);

   real dampk3 = dampk * dampk2;
   dmpk[0] = 1 - (1 + 0.5f * dampk) * expk;
   dmpk[1] = 1 - (1 + dampk + 0.5f * dampk2) * expk;
   dmpk[2] = 1 - (1 + dampk + 0.5f * dampk2 + dampk3 * div6) * expk;
   //} // endif GORDON1


   // if CONSTEXPR (eq<PENTYP, GORDON2>()) {
   //    const real div3 = 1 / ((real)3);

   //    dmpk[0] = 1 - expk;
   //    dmpk[1] = 1 - (1 + dampk) * expk;
   //    dmpk[2] = 1 - (1 + dampk + dampk2 * div3) * expk;
   // } // endif GORDON2
}
}
