#pragma once
#include "mathfunc.h"
#include "named_struct.h"
#include "seq_def.h"


TINKER_NAMESPACE_BEGIN
#pragma acc routine seq
template <class Ver, class ETYP>
SEQ_CUDA
void pair_charge(real r, real xr, real yr, real zr, real cscale, real chgi,
                 real chgk, real ebuffer, real f, real aewald,
                 real& restrict grdx, real& restrict grdy, real& restrict grdz,
                 int& restrict ctl, real& restrict etl, real& restrict vtlxx,
                 real& restrict vtlxy, real& restrict vtlxz,
                 real& restrict vtlyy, real& restrict vtlyz,
                 real& restrict vtlzz)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;
   MAYBE_UNUSED real e, dedx, dedy, dedz;


   if CONSTEXPR (eq<ETYP, EWALD>()) {
      real fik = f * chgi * chgk;
      real rew = aewald * r;
      real erfterm = REAL_ERFC(rew);
      real rb = r + ebuffer;
      real invrb = REAL_RECIP(rb);


      if CONSTEXPR (do_e) {
         e = fik * invrb * erfterm;
      }
      if CONSTEXPR (do_g) {
         real invr = REAL_RECIP(r);
         real invrb2 = invrb * invrb;
         real de = erfterm * invrb2;
         de += (2 * aewald / sqrtpi) * REAL_EXP(-rew * rew) * invr;
         de *= -fik;
         de *= invr;
         dedx = de * xr;
         dedy = de * yr;
         dedz = de * zr;
      }
   } else if CONSTEXPR (eq<ETYP, NON_EWALD>()) {
      real fik = cscale * f * chgi * chgk;
      real rb = r + ebuffer;
      real invrb = REAL_RECIP(rb);


      if CONSTEXPR (do_e) {
         e = fik * invrb;
      }
      if CONSTEXPR (do_g) {
         real invr = REAL_RECIP(r);
         real invrb2 = invrb * invrb;
         real de = -fik * invrb2;
         de *= invr;
         dedx = de * xr;
         dedy = de * yr;
         dedz = de * zr;
      }
   }


   if CONSTEXPR (do_a) {
      // ci and ck may be zero
      if (e != 0) {
         // scale =  1.0 -- count +1  scale =  0.4 or 0.0 -- impossible
         // scale = -0.6 -- count +0  scale = -1.0        -- count -1
         if (cscale == 1)
            ctl += 1;
         else if (cscale == -1)
            ctl -= 1;
      }
   }
   if CONSTEXPR (do_e) {
      etl += e;
   }
   if CONSTEXPR (do_g) {
      grdx += dedx;
      grdy += dedy;
      grdz += dedz;
   }
   if CONSTEXPR (do_v) {
      real vxx = xr * dedx;
      real vyx = yr * dedx;
      real vzx = zr * dedx;
      real vyy = yr * dedy;
      real vzy = zr * dedy;
      real vzz = zr * dedz;
      vtlxx += vxx;
      vtlxy += vyx;
      vtlxz += vzx;
      vtlyy += vyy;
      vtlyz += vzy;
      vtlzz += vzz;
   }
}
TINKER_NAMESPACE_END
