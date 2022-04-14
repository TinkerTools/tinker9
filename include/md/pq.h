#pragma once
#include "ff/atom.h"

namespace tinker {
/// \ingroup mdpq
/// \brief Update #x, #y, #z via `x += v * dt`.
/// Currently #xpos, #ypos, and #zpos are integrated first, then uses
/// #copyPosToXyz() to update #x, #y, #z.
void mdPos(time_prec dt);
void mdPos(time_prec, pos_prec*, pos_prec*, pos_prec*, //
   const vel_prec*, const vel_prec*, const vel_prec*);

/// \ingroup mdpq
/// \brief Update #xp, #y, #z via `x = a x + b v`.
void mdPosAxbv(pos_prec a, pos_prec b);

/// \ingroup mdpq
/// \brief Update #xp, #y, #z via `r = MatA R + MatB v`.
void mdPosAxbvAn(pos_prec (*a)[3], pos_prec (*b)[3]);

/// \ingroup mdpq
void mdDebugPosNorm(pos_prec poseps, time_prec dt, //
   const vel_prec* vx, const vel_prec* vy, const vel_prec* vz);

/// \ingroup mdpq
/// \brief Update velocities via `v += -g/m dt`.
void mdVel(time_prec dt, const grad_prec* grx, const grad_prec* gry, const grad_prec* grz);

/// \ingroup mdpq
/// \brief Update velocities via `v = v0 -g/m dt`.
void mdVelB(time_prec dt, vel_prec* vlx, vel_prec* vly, vel_prec* vlz, //
   const vel_prec* vlx0, const vel_prec* vly0, const vel_prec* vlz0,   //
   const grad_prec* grx, const grad_prec* gry, const grad_prec* grz);

/// \ingroup mdpq
/// \brief Update velocities via `v += (-g/m dt -g2/m dt2)`.
void mdVel2(time_prec dt, const grad_prec* grx, const grad_prec* gry, const grad_prec* grz,
   time_prec dt2, const grad_prec* grx2, const grad_prec* gry2, const grad_prec* grz2);

/// \ingroup mdpq
void mdVelScale(vel_prec scal, int nelem, vel_prec* vx, vel_prec* vy, vel_prec* vz);

/// \ingroup mdpq
/// \brief Update velocities via `v = a v + b (g1/nrespa + g2)/m t` (isotropic).
/// If \p nrespa equals 1, the 2nd set of gradients are ignored.
void mdVelAvbf(int nrespa, vel_prec a, vel_prec b,                   //
   const grad_prec* gx1, const grad_prec* gy1, const grad_prec* gz1, //
   const grad_prec* gx2, const grad_prec* gy2, const grad_prec* gz2);

/// \ingroup mdpq
/// \brief Update velocities via `v = matA v + matB (g1/nrespa + g2)/m t` (anisotropic).
/// If \p nrespa equals 1, the 2nd set of gradients are ignored.
/// \see #mdVelAvbf.
void mdVelAvbfAn(int nrespa, vel_prec a[3][3], vel_prec b[3][3],     //
   const grad_prec* gx1, const grad_prec* gy1, const grad_prec* gz1, //
   const grad_prec* gx2, const grad_prec* gy2, const grad_prec* gz2);

/// \ingroup mdpq
void mdVelData(RcOp);
}

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

namespace tinker {
/// \ingroup mdpq
/// \{
/// \var vx
/// \brief Velocities.
/// \var vy
/// \copydoc vx
/// \var vz
/// \copydoc vx
/// \}
TINKER_EXTERN vel_prec *vx, *vy, *vz;
}
