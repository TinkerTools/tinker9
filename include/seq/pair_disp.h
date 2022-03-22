#pragma once
#include "math/inc.h"
#include "seq/damp_hippodisp.h"
#include "seq/switch.h"
#include "seqdef.h"

namespace tinker {
#pragma acc routine seq
template <bool DO_G, class DTYP, int SCALE>
SEQ_CUDA
void pair_disp_obsolete(real r, real r2, real rr1, //
   real dspscale, real aewald, real ci, real ai, real ck, real ak, real edcut, real edoff,
   real& restrict e, real& restrict de)
{
   if (r > edoff) {
      e = 0;
      if CONSTEXPR (DO_G) {
         de = 0;
      }
      return;
   }

#if TINKER_REAL_SIZE == 8
   real eps = 0.001f;
#elif TINKER_REAL_SIZE == 4
   real eps = 0.05f;
#endif

   real diff = REAL_ABS(ai - ak);
   real rr2 = rr1 * rr1;
   real rr6 = rr2 * rr2 * rr2;
   real di = ai * r;
   real dk = ak * r;
   real expi = REAL_EXP(-di);
   real expk = REAL_EXP(-dk);
   real di2 = di * di;
   real damp = 1, ddamp = 0;
   if (diff > eps) {
      real ai2 = ai * ai;
      real ak2 = ak * ak;
      real ti = ak2 * REAL_RECIP((ak + ai) * (ak - ai));
      real tk = ai2 * REAL_RECIP((ai + ak) * (ai - ak));
      real a1 = 2 * ti + 1;
      real b1 = 2 * tk + 1;
      real termi = ((di / 2 + b1) / 2 * di + b1) * di + b1;
      real termk = ((dk / 2 + a1) / 2 * dk + a1) * dk + a1;
      termi *= ti * ti * expi;
      termk *= tk * tk * expk;
      damp = 1 - termi - termk;
      if CONSTEXPR (DO_G) {
         real dk2 = dk * dk;
         real a2 = (di - 1) / 4 + tk;
         real b2 = (dk - 1) / 4 + ti;
         a2 *= ai * di2 * ti * ti * expi;
         b2 *= ak * dk2 * tk * tk * expk;
         ddamp = a2 + b2;
      }
   } else {
      ai = (ai + ak) * 0.5f;
      di = ai * r;
      di2 = di * di;
      expi = REAL_EXP(-di);
      real term = ((((di + 5) * di + 17) * di / 96 + 0.5f) * di + 1) * di + 1;
      damp = 1 - term * expi;
      if CONSTEXPR (DO_G)
         ddamp = ai * expi * di2 * ((di2 - 3) * di - 3) / 96;
   }

   if CONSTEXPR (SCALE == 1)
      dspscale = 1;

   if CONSTEXPR (eq<DTYP, DEWALD>()) {
      real ralpha2 = r2 * aewald * aewald;
      real term = 1 + ralpha2 + 0.5f * ralpha2 * ralpha2;
      real expterm = REAL_EXP(-ralpha2);
      real expa = expterm * term;
      e = -ci * ck * rr6 * (dspscale * damp * damp + expa - 1);
      if CONSTEXPR (DO_G) {
         real rterm = -ralpha2 * ralpha2 * ralpha2 * rr1 * expterm;
         de = -6 * e * rr1 - ci * ck * rr6 * (rterm + 2 * dspscale * damp * ddamp);
      }
   } else if CONSTEXPR (eq<DTYP, NON_EWALD_TAPER>()) {
      e = -ci * ck * rr6;
      if CONSTEXPR (DO_G) {
         de = -6 * e * rr1;
         de = de * damp * damp + 2 * e * damp * ddamp;
      }
      e = e * damp * damp;
      if (r > edcut) {
         real taper, dtaper;
         switch_taper5<DO_G>(r, edcut, edoff, taper, dtaper);
         if CONSTEXPR (DO_G)
            de = e * dtaper + de * taper;
         e = e * taper;
      }
      e *= dspscale;
      if CONSTEXPR (DO_G)
         de *= dspscale;
   }
}

#pragma acc routine seq
template <bool DO_G, class DTYP, int SCALE>
SEQ_CUDA
void pair_disp(real r, real r2, real rr1, //
   real dspscale, real aewald, real ci, real ai, real ck, real ak, real edcut, real edoff,
   real& restrict e, real& restrict de)
{
   if (r > edoff) {
      e = 0;
      if CONSTEXPR (DO_G) {
         de = 0;
      }
      return;
   }

   real rr2 = rr1 * rr1;
   real rr6 = rr2 * rr2 * rr2;
   real dmpik[2], damp, ddamp;
   damp_hippodisp<DO_G>(dmpik, r, rr1, ai, ak);
   damp = dmpik[0];
   if CONSTEXPR (DO_G)
      ddamp = dmpik[1];

   if CONSTEXPR (SCALE == 1)
      dspscale = 1;
   if CONSTEXPR (eq<DTYP, DEWALD>()) {
      real ralpha2 = r2 * aewald * aewald;
      real term = 1 + ralpha2 + 0.5f * ralpha2 * ralpha2;
      real expterm = REAL_EXP(-ralpha2);
      real expa = expterm * term;
      e = -ci * ck * rr6 * (dspscale * damp * damp + expa - 1);
      if CONSTEXPR (DO_G) {
         real rterm = -ralpha2 * ralpha2 * ralpha2 * rr1 * expterm;
         de = -6 * e * rr1 - ci * ck * rr6 * (rterm + 2 * dspscale * damp * ddamp);
      }
   } else if CONSTEXPR (eq<DTYP, NON_EWALD_TAPER>()) {
      e = -ci * ck * rr6;
      if CONSTEXPR (DO_G) {
         de = -6 * e * rr1;
         de = de * damp * damp + 2 * e * damp * ddamp;
      }
      e = e * damp * damp;
      if (r > edcut) {
         real taper, dtaper;
         switch_taper5<DO_G>(r, edcut, edoff, taper, dtaper);
         if CONSTEXPR (DO_G)
            de = e * dtaper + de * taper;
         e = e * taper;
      }
      e *= dspscale;
      if CONSTEXPR (DO_G)
         de *= dspscale;
   }
}
}
