#pragma once
#include "ff/evdw.h"
#include "seq/seq.h"

namespace tinker {
/**
 * \ingroup vdw
 */
SEQ_ROUTINE
inline real pair_vlambda(real vlam, Vdw vcouple, int imut, int kmut)
{
   real vlambda = 1;
   if (vcouple == Vdw::DECOUPLE) {
      vlambda = (imut == kmut ? 1 : vlam);
   } else if (vcouple == Vdw::ANNIHILATE) {
      vlambda = (imut || kmut ? vlam : 1);
   }
   return vlambda;
}
}
