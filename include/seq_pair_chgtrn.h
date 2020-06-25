#pragma once
#include "mathfunc.h"
#include "seq_def.h"


namespace tinker {
#pragma acc routine seq
template <bool DO_G>
SEQ_CUDA
void pair_chgtrn(real r, real mscale, real f, real alphai, real chgi,
                 real alphak, real chgk, e_prec& restrict e,
                 e_prec& restrict de)
{
   f *= mscale;
   real expi = REAL_EXP(-alphai * r);
   real expk = REAL_EXP(-alphak * r);
   e = -chgi * expk - chgk * expi;
   e *= f;
   if CONSTEXPR (DO_G) {
      de = chgi * expk * alphak + chgk * expi * alphai;
      de *= f;
   }
}
}
