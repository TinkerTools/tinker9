#pragma once
#include "compare_types.h"


namespace tinker {
namespace calc {
/**
 * \ingroup mdcalc
 * \brief Bitmasks for MD.
 */
enum __CalcBits : int
{
   xyz = 0x001,  ///< Use coordinates.
   vel = 0x002,  ///< Use velocities.
   mass = 0x004, ///< Use mass.
   traj = 0x008, ///< Use multi-frame trajectory.


   energy = 0x010, ///< Evaluate energy.
   grad = 0x020,   ///< Evaluate energy gradient.
   virial = 0x040, ///< Evaluate virial tensor.
   analyz = 0x080, ///< Evaluate number of interactions.


   vmask = energy + grad + virial +
      analyz,   ///< Bits mask to clear energy-irrelevant flags.
   v0 = energy, ///< Similar to Tinker energy routines. Energy only.
   v1 = energy + grad + virial, ///< Similar to version 1 Tinker energy
                                ///< routines. Energy, gradient, and virial.
   v3 = energy + analyz, ///< Similar to version 3 Tinker energy routines.
                         ///< Energy and number of interactions.
   v4 = energy + grad,   ///< Energy and gradient.
   v5 = grad,            ///< Gradient only.
   v6 = grad + virial,   ///< Gradient and virial.


   md = 0x100, ///< Run MD simulation.
};


/**
 * \ingroup mdcalc
 * \brief Sanity checks for version constants.
 */
template <int USE>
struct __Vers
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


extern "C"
{
   // PME grids.
   struct PCHG
   {};
   struct MPOLE
   {};
   struct UIND
   {};
   struct UIND2
   {};
   struct DISP
   {};


   // Ewald vs. Non-Ewald
   struct EWALD
   {};
   struct DEWALD
   {}; // Dispersion PME
   struct NON_EWALD
   {};
   struct NON_EWALD_TAPER
   {}; // Non-EWALD partial charge also uses switching functions.


   // Energy versions.
   struct Eng : public tinker::calc::__Vers<tinker::calc::v0>
   {};
   struct EngGradVir : public tinker::calc::__Vers<tinker::calc::v1>
   {};
   struct EngAlyz : public tinker::calc::__Vers<tinker::calc::v3>
   {};
   struct EngGrad : public tinker::calc::__Vers<tinker::calc::v4>
   {};
   struct Grad : public tinker::calc::__Vers<tinker::calc::v5>
   {};
   struct GradVir : public tinker::calc::__Vers<tinker::calc::v6>
   {};


   // RATTLE vs SHAKE
   struct RATTLE
   {};
   struct SHAKE
   {};


   // Bond terms.
   struct HARMONIC
   {};
   struct MORSE
   {};


   // Opbend terms.
   struct WDC
   {};
   struct ALLINGER
   {};


   // VDW terms.
   struct LJ
   {};
   struct BUCK
   {};
   struct MM3HB
   {};
   struct HAL
   {};
   struct GAUSS
   {};
}


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
