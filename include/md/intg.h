#pragma once
#include "ff/precision.h"
#include "ff/timescale.h"
#include "tool/rcman.h"

namespace tinker {
void mdData(RcOp);
void mdIntegrateData(RcOp);

void mdrest(int istep);
void mdPropagate(int nsteps, time_prec dt_ps);

constexpr unsigned RESPA_FAST = 1; // 2**0, fast group shall be 0.
constexpr unsigned RESPA_SLOW = 2; // 2**1, slow group shall be 1.
const TimeScaleConfig& respaTSConfig();

void mdsaveAsync(int istep, time_prec dt);
void mdsaveSynchronize();
void mdsaveData(RcOp);
}

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

namespace tinker {
/// \ingroup mdintg
/// \{
/// \var gx1
/// \brief Gradient for the fast RESPA energy terms.
/// \var gy1
/// \copydoc gx1
/// \var gz1
/// \copydoc gx1
/// \var gx2
/// \brief Gradient for the slow RESPA energy terms.
/// \var gy2
/// \copydoc gx2
/// \var gz2
/// \copydoc gx2
/// \}
TINKER_EXTERN grad_prec *gx1, *gy1, *gz1;
TINKER_EXTERN grad_prec *gx2, *gy2, *gz2;
}
