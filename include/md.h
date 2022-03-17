#pragma once
#include "precision.h"

//====================================================================//
// mdcalc

#include "compare_types.h"

namespace tinker {
namespace calc {
/// \ingroup mdcalc
/// \brief Bitmasks for MD.
enum Bits : int
{
   xyz = 0x001,  ///< Use coordinates.
   vel = 0x002,  ///< Use velocities.
   mass = 0x004, ///< Use mass.
   traj = 0x008, ///< Use multi-frame trajectory.

   energy = 0x010, ///< Evaluate energy.
   grad = 0x020,   ///< Evaluate energy gradient.
   virial = 0x040, ///< Evaluate virial tensor.
   analyz = 0x080, ///< Evaluate number of interactions.
   md = 0x100,     ///< Run MD simulation.

   vmask = energy + grad + virial + analyz, ///< Bits mask to clear energy-irrelevant flags.
   v0 = energy,                             ///< Similar to Tinker energy routines. Energy only.
   v1 = energy + grad + virial,             ///< Similar to version 1 Tinker energy
                                            ///< routines. Energy, gradient, and virial.
   v3 = energy + analyz,                    ///< Similar to version 3 Tinker energy routines.
                                            ///< Energy and number of interactions.
   v4 = energy + grad,                      ///< Energy and gradient.
   v5 = grad,                               ///< Gradient only.
   v6 = grad + virial,                      ///< Gradient and virial.
};

/// \ingroup mdcalc
/// \brief Sanity checks for version constants.
template <int USE>
struct Vers
{
   static constexpr int value = USE;
   static constexpr int e = USE & calc::energy;
   static constexpr int a = USE & calc::analyz;
   static constexpr int g = USE & calc::grad;
   static constexpr int v = USE & calc::virial;
   static_assert(v ? (bool)g : true, "If calc::virial, must calc::grad.");
   static_assert(a ? (bool)e : true, "If calc::analyz, must calc::energy.");
};
}
}

#include "cname.h"

namespace tinker {
namespace calc {
using V0 = Eng;
using V1 = EngGradVir;
using V3 = EngAlyz;
using V4 = EngGrad;
using V5 = Grad;
using V6 = GradVir;
}
}

//====================================================================//
// mdpq

#include "tool/rc_man.h"
#include <istream>

namespace tinker {
/// \ingroup mdpq
void mdNData(rc_op);

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
void mdXyzData(rc_op);

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
void mdVelData(rc_op);

/// \ingroup mdpq
void mdMassData(rc_op);
}

//====================================================================//
// mdpt

namespace tinker {
void mdKineticEnergy(energy_prec& eksum_out, energy_prec (&ekin_out)[3][3], int n,
   const double* mass, const vel_prec* vx, const vel_prec* vy, const vel_prec* vz);
void mdKineticExplicit(T_prec& temp_out, energy_prec& eksum_out, energy_prec (&ekin_out)[3][3],
   const vel_prec* vx, const vel_prec* vy, const vel_prec* vz);
void mdKinetic(T_prec& temp);

void mdBussiThermostat(time_prec dt, T_prec temp);
void mdBussiThermostat_acc(time_prec dt, T_prec temp);

/// \ingroup mdpt
/// \brief Applies a box size correction as needed for the Monte Carlo barostat
/// at the half time step.
///
/// Literature reference:
///    - <a href="https://doi.org/10.1080/00268977200100031">
///    I. R. McDonald,
///    "NpT-ensemble Monte Carlo calculations for binary liquid mixtures",
///    Molecular Physics, 23, 41-58 (1972).
///    </a>
void mdMonteCarloBarostat(energy_prec epot, T_prec temp);
void mdMonteCarloBarostat_acc(energy_prec epot, T_prec temp);

/// \ingroup mdpt
/// \brief Berendsen barostat by scaling the coordinates and box dimensions via
/// coupling to an external constant pressure bath. Code for anisotropic pressure
/// coupling was provided by Guido Raos, Dipartimento di Chimica, Politecnico di
/// Milano, Italy.
///
/// Literature reference:
///    - <a href="https://doi.org/10.1063/1.448118">
///    H. J. C. Berendsen, J. P. M. Postma, W. F. van Gunsteren, A. DiNola,
///    and J. R. Hauk,
///    "Molecular dynamics with coupling to an external bath",
///    J. Chem. Phys., 81, 3684-3690 (1984).
///    </a>
void mdBerendsenBarostat(time_prec dt);
void mdBerendsenBarostat_acc(time_prec);
}

//====================================================================//
// mdintg

#include "time_scale.h"

namespace tinker {
void mdrest(int istep);
void mdData(rc_op op);
void mdPropagate(int nsteps, time_prec dt_ps);
void mdIntegrateData(rc_op);

constexpr unsigned RESPA_FAST = 1; // 2**0, fast group shall be 0.
constexpr unsigned RESPA_SLOW = 2; // 2**1, slow group shall be 1.
const TimeScaleConfig& mdRespaTsconfig();

void mdsaveAsync(int istep, time_prec dt);
void mdsaveSynchronize();
void mdsaveData(rc_op);
}

#include "egv.h"

#include "glob.md.h"
#include "glob.nelembuffer.h"
