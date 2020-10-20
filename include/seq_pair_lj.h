#pragma once
#include "mathfunc.h"
#include "seq_def.h"
#include "seq_switch.h"


namespace tinker {
/**
 * \ingroup vdw
 */
#pragma acc routine seq
template <bool DO_G>
SEQ_CUDA
void pair_lj(real rik, real rik2, real rv, real eps, real vscalek,
             real& restrict e, real& restrict de)
{
   eps *= vscalek;
   real rv2_rik2 = rv * rv * REAL_RECIP(rik2);
   real p6 = rv2_rik2 * rv2_rik2 * rv2_rik2;
   e = eps * p6 * (p6 - 2);
   if CONSTEXPR (DO_G) {
      de = eps * p6 * (p6 - 1) * (-12 * REAL_RECIP(rik));
   }
}


/**
 * \ingroup vdw
 */
#pragma acc routine seq
template <bool DO_G, class RADRULE, class EPSRULE, int SCALE>
SEQ_CUDA
void pair_lj_v2(real r, real invr, //
                real vscale, real radi, real epsi, real radk, real epsk,
                real evcut, real evoff, real& restrict ev, real& restrict dev)
{
   if (r > evoff) {
      ev = 0;
      if CONSTEXPR (DO_G) {
         dev = 0;
      }
      return;
   }
   real rad = RADRULE::avg(radi, radk);
   real eps = EPSRULE::savg(epsi, epsk);
   if CONSTEXPR (SCALE != 1)
      eps *= vscale;
   real sig = invr * rad;
   real sig2 = sig * sig;
   real p6 = sig2 * sig2 * sig2;
   real ep6 = eps * p6;
   ev = ep6 * (p6 - 2);
   if CONSTEXPR (DO_G) {
      dev = ep6 * (p6 - 1) * (-12 * invr);
   }
   // taper
   if (r > evcut) {
      real taper, dtaper;
      switch_taper5<DO_G>(r, evcut, evoff, taper, dtaper);
      if CONSTEXPR (DO_G)
         dev = ev * dtaper + dev * taper;
      ev = ev * taper;
   }
}


/**
 * \ingroup vdw
 */
#pragma acc routine seq
template <bool DO_G, int SCALE>
SEQ_CUDA
void pair_lj_v3(real r, real invr, //
                real vscale, real rad, real eps, real evcut, real evoff,
                real& restrict ev, real& restrict dev)
{
   if CONSTEXPR (SCALE != 1) {
      eps *= vscale;
   }
   real sig = invr * rad;
   real sig2 = sig * sig;
   real p6 = sig2 * sig2 * sig2;
   real ep6 = eps * p6;
   ev = ep6 * (p6 - 2);
   if CONSTEXPR (DO_G) {
      dev = ep6 * (p6 - 1) * (-12 * invr);
   }
   // taper
   if (r > evcut) {
      real taper, dtaper;
      switch_taper5<DO_G>(r, evcut, evoff, taper, dtaper);
      if CONSTEXPR (DO_G)
         dev = ev * dtaper + dev * taper;
      ev = ev * taper;
   }
}
}
