#pragma once
#include "ff/atom.h"
#include <istream>

namespace tinker {
/// \ingroup mdpq
/// \brief Update #x, #y, #z by #xpos, #ypos, and #zpos.
/// If #xpos etc. are only aliases, return directly.
void mdCopyPosToXyz();
void mdCopyPosToXyz(bool check_nblist);
void mdCopyPosToXyz_acc();

/// \ingroup mdpq
/// \brief Update #x, #y, #z via `x += v * dt`.
/// Currently #xpos, #ypos, and #zpos are integrated first, then uses
/// #mdCopyPosToXyz() to update #x, #y, #z.
void mdPos(time_prec dt);
void mdPos(time_prec, pos_prec*, pos_prec*, pos_prec*, //
   const vel_prec*, const vel_prec*, const vel_prec*);
void mdPos_acc(time_prec, pos_prec*, pos_prec*, pos_prec*, //
   const vel_prec*, const vel_prec*, const vel_prec*);

/// \ingroup mdpq
/// \brief Update #xp, #y, #z via `x = a x + b v`.
void mdPosAxbv(pos_prec a, pos_prec b);

/// \ingroup mdpq
void mdDebugPosNorm_acc(pos_prec poseps, time_prec dt, //
   const vel_prec* vx, const vel_prec* vy, const vel_prec* vz);

/// \ingroup mdpq
/// \brief Call #mdBounds() at least every x steps in MD.
constexpr int BOUNDS_EVERY_X_STEPS = 500;

/// \ingroup mdpq
/// \brief Finds the geometric center of each molecule and translate any stray
/// molecules back into the periodic box on GPU.
/// \note
///    - Updating #x, #y, #z is the goal.
///    - Checks whether PBC is in use inside this function.
///    - Will not perturb the neighbor lists so no need to update them.
///    - Tinker uses centers of mass.
void mdBounds();
void mdBounds_acc();

/// \ingroup mdpq
void mdReadFrameCopyinToXyz(std::istream& input, int& done);

/// \ingroup mdpq
/// \brief Update velocities via `v += -g/m dt`.
void mdVel(time_prec dt, const grad_prec* grx, const grad_prec* gry, const grad_prec* grz);
void mdVelA_acc(time_prec, vel_prec*, vel_prec*, vel_prec*, //
   const grad_prec*, const grad_prec*, const grad_prec*);

/// \ingroup mdpq
/// \brief Update velocities via `v = v0 -g/m dt`.
void mdVelB(time_prec dt, vel_prec* vlx, vel_prec* vly, vel_prec* vlz, //
   const vel_prec* vlx0, const vel_prec* vly0, const vel_prec* vlz0,   //
   const grad_prec* grx, const grad_prec* gry, const grad_prec* grz);
void mdVelB_acc(time_prec, vel_prec*, vel_prec*, vel_prec*, //
   const vel_prec*, const vel_prec*, const vel_prec*,       //
   const grad_prec*, const grad_prec*, const grad_prec*);

/// \ingroup mdpq
/// \brief Update velocities via `v += (-g/m dt -g2/m dt2)`.
void mdVel2(time_prec dt, const grad_prec* grx, const grad_prec* gry, const grad_prec* grz,
   time_prec dt2, const grad_prec* grx2, const grad_prec* gry2, const grad_prec* grz2);
void mdVel2_acc(time_prec, const grad_prec*, const grad_prec*, const grad_prec*, time_prec,
   const grad_prec*, const grad_prec*, const grad_prec*);

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
