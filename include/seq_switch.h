#pragma once
#include "macro_void_cuda_def.h"
#include "mathfunc.h"
#include "switch.h"


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup ff
 * \brief Smooth function
 * `F: [cut,off]->[1,0]`.
 *
 * Derived from the second order smooth step function
 * \f[ S_2: [0,1]\rightarrow[0,1] \f]
 * \f[ S_2(x) = 6 x^5 - 15 x^4 + 10 x^3 \f]
 *
 * \param[in] rik     Distance.
 * \param[in] cut     Distance at which switching of the potential begins.
 * \param[out] off    Distance at which the potential energy goes to zero.
 * \param[out] taper  F value.
 * \param[out] dtaper dF/dx value.
 * \tparam DO_DTAPER  If false, `dtaper` will not be calculated and its
 *                    original value will not be altered.
 */
#pragma acc routine seq
template <int DO_DTAPER>
__device__
void switch_taper5(real rik, real cut, real off, real& restrict taper,
                   real& restrict dtaper)
{
   real rinv = REAL_RECIP(cut - off);
   real x = (rik - off) * rinv;
   real x2 = x * x;
   real x3 = x2 * x;
   taper = x3 * (6 * x2 - 15 * x + 10);
   if CONSTEXPR (DO_DTAPER) {
      dtaper = 30 * (x * (1 - x)) * (x * (1 - x)) * rinv;
   }
}
TINKER_NAMESPACE_END
