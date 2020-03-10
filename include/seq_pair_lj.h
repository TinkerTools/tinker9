#pragma once
#include "mathfunc.h"
#include "seq_def.h"


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup vdw
 */
#pragma acc routine seq
template <int DO_G>
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
TINKER_NAMESPACE_END
