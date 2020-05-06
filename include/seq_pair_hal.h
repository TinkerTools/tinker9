#pragma once
#include "mathfunc.h"
#include "seq_def.h"


namespace tinker {
/**
 * \ingroup vdw
 */
#pragma acc routine seq
template <int DO_G>
SEQ_CUDA
void pair_hal(real rik, real rv, real eps, real vscalek, real vlambda, //
              real ghal, real dhal, real scexp, real scalpha,          //
              real& restrict e, real& restrict de)
{
   eps *= vscalek;
   real rho = rik * REAL_RECIP(rv);
   real rho6 = REAL_POW(rho, 6);
   real rho7 = rho6 * rho;
   eps *= REAL_POW(vlambda, scexp);
   real scal = scalpha * REAL_POW(1 - vlambda, 2);
   real s1 = REAL_RECIP(scal + REAL_POW(rho + dhal, 7));
   real s2 = REAL_RECIP(scal + rho7 + ghal);
   real t1 = REAL_POW(1 + dhal, 7) * s1;
   real t2 = (1 + ghal) * s2;
   e = eps * t1 * (t2 - 2);
   if CONSTEXPR (DO_G) {
      real dt1drho = -7 * REAL_POW(rho + dhal, 6) * t1 * s1;
      real dt2drho = -7 * rho6 * t2 * s2;
      de = eps * (dt1drho * (t2 - 2) + t1 * dt2drho) * REAL_RECIP(rv);
   }
}
}
