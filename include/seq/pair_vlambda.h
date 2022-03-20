#pragma once
#include "glob/vdwpot.h"
#include "seq/switch.h"
#include "seq_def.h"

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
