#pragma once
#include "ff/pchg/evdw1.h"
#include "math/switch.h"
#include "seqdef.h"

namespace tinker {
/**
 * \ingroup vdw
 */
SEQ_ROUTINE
inline real pair_vlambda(real vlam, evdw_t vcouple, int imut, int kmut)
{
   real vlambda = 1;
   if (vcouple == evdw_t::decouple) {
      vlambda = (imut == kmut ? 1 : vlam);
   } else if (vcouple == evdw_t::annihilate) {
      vlambda = (imut || kmut ? vlam : 1);
   }
   return vlambda;
}
}
