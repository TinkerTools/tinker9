#pragma once
#include "mathfunc.h"
#include "seq_def.h"


namespace tinker {
#pragma acc routine seq
template <bool DO_G, class DTYP>
SEQ_CUDA
void pair_disp(real r, real r2, real dspscale, real aewald, real ci, real ai,
               real ck, real ak, e_prec& restrict e, e_prec& restrict de)
{
   real rr1 = REAL_RECIP(r);
   real rr2 = REAL_RECIP(r2);
   real rr6 = rr2 * rr2 * rr2;


   real di = ai * r;
   real dk = ak * r;
   real expi = REAL_EXP(-di);
   real expk = REAL_EXP(-dk);
   real di2 = di * di;
   real damp = 1, ddamp = 0;
   if (ai != ak) {
      real ai2 = ai * ai;
      real ak2 = ak * ak;
      real ti = ak2 * REAL_RECIP(ak2 - ai2);
      real tk = ai2 * REAL_RECIP(ai2 - ak2);
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
   }
   if (ai == ak) {
      real term = ((((di + 5) * di + 17) * di / 96 + 0.5f) * di + 1) * di + 1;
      damp = 1 - term * expi;
      if CONSTEXPR (DO_G)
         ddamp = ai * expi * di2 * ((di2 - 3) * di - 3) / 96;
   }


   if CONSTEXPR (eq<DTYP, DEWALD>()) {
      real ralpha2 = r2 * aewald * aewald;
      real term = 1 + ralpha2 + 0.5f * ralpha2 * ralpha2;
      real expterm = REAL_EXP(-ralpha2);
      real expa = expterm * term;
      e = -ci * ck * rr6 * (expa + damp * damp - 1);
      if CONSTEXPR (DO_G) {
         real rterm = -ralpha2 * ralpha2 * ralpha2 * rr1 * expterm;
         de = -6 * e * rr1 - ci * ck * rr6 * (rterm + 2 * damp * ddamp);
      }
   } else if CONSTEXPR (eq<DTYP, NON_EWALD>() || eq<DTYP, NON_EWALD_TAPER>()) {
      e = -ci * ck * rr6;
      if CONSTEXPR (DO_G) {
         de = -6 * e * rr1;
         de = dspscale * (de * damp * damp + 2 * e * damp * ddamp);
      }
      e = dspscale * (e * damp * damp);
   }
}
}
