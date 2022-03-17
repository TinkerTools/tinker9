#pragma once
#include "macro.h"

//====================================================================//
// mdprec

/**
 * \ingroup prec
 * \page prec
 *
 * | Properties                  | Types       | Underlying Types |
 * |----------------------------:|------------:|-----------------:|
 * | time                        | time_prec   | mixed            |
 * | temperature                 | T_prec      | mixed            |
 * | velocity                    | vel_prec    | mixed            |
 * | position (for integrators)  | pos_prec    | mixed            |
 * | energy (for integrators)    | energy_prec | mixed            |
 * | virial (for integrators)    | virial_prec | mixed            |
 * | gradient (deterministic)    | grad_prec   | fixed            |
 * | gradient (floating-point)   | grad_prec   | real             |
 * | energy (individual terms)   | e_prec      | real             |
 * | virial (individual terms)   | v_prec      | real             |
 * | gradient (individual terms) | real        | real             |
 * | position (for energies)     | real        | real             |
 */

/**
 * \ingroup prec
 * \def TINKER_DETERMINISTIC_FORCE
 * \brief Logical macro for the underlying type of energy gradients.
 *    - If `true`, always use fixed-point arithmetic to accumulate energy
 *    gradients, regardless of the underlying type of each individual component.
 *    - If `false`, use type `real`, which can either be `float` or `double`
 *    based on the precision macro.
 *
 * \see TINKER_DOUBLE_PRECISION
 * \see TINKER_MIXED_PRECISION
 * \see TINKER_SINGLE_PRECISION
 *
 * In general, evaluating energy, forces, etc. twice, we don't expect to get
 * two identical answers, but we may not care as much because the difference
 * is usually negligible.
 * [See [Why is cos(x) !=
 * cos(y)?](https://isocpp.org/wiki/faq/newbie#floating-point-arith2)]
 * Whereas in MD, two simulations with the same initial configurations can
 * easily diverge due to the accumulated difference. If, for whatever reason,
 * you are willing to elongate the process of the inevitable divergence at the
 * cost of slightly slower simulation speed, a more "deterministic" force (using
 * fixed-point arithmetic) can help.
 *
 * To accumulate floating-point values via fixed-point arithmetic, we first
 * scale the values by a large integer then sum only the integer part.
 * \code{.cpp}
 * // fixed sum;
 * // FLT is float or double
 * FLT val2 = val * 0x100000000ull;
 * sum += (fixed)((long long)val2);
 * \endcode
 * To get the floating-point sum, we have to cast the integer sum back to
 * floating-point number, then divide it by the same large integer.
 * \code{.cpp}
 * FLT val2 = (FLT)((long long)sum);
 * FLT answer = val2 / 0x100000000ull;
 * \endcode
 */
#ifndef TINKER_DETERMINISTIC_FORCE
#   if TINKER_DOUBLE_PRECISION
#      define TINKER_DETERMINISTIC_FORCE 0
#   else
#      define TINKER_DETERMINISTIC_FORCE 1
#   endif
#endif

namespace tinker {
using time_prec = mixed;
using T_prec = mixed;
using vel_prec = mixed;
using pos_prec = mixed;
using energy_prec = mixed; // total energies
using virial_prec = mixed; // total virial tensor
using e_prec = real;       // individual energy
using v_prec = real;       // individual virial
#if TINKER_DETERMINISTIC_FORCE
using grad_prec = fixed;
#else
using grad_prec = real;
#endif
}

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
   struct Eng : public tinker::calc::Vers<tinker::calc::v0>
   {};
   struct EngGradVir : public tinker::calc::Vers<tinker::calc::v1>
   {};
   struct EngAlyz : public tinker::calc::Vers<tinker::calc::v3>
   {};
   struct EngGrad : public tinker::calc::Vers<tinker::calc::v4>
   {};
   struct Grad : public tinker::calc::Vers<tinker::calc::v5>
   {};
   struct GradVir : public tinker::calc::Vers<tinker::calc::v6>
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

   // GORDON1 vs. GORDON2 damping functions
   struct GORDON1
   {};
   struct GORDON2
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

//====================================================================//
// mdpq
// mdegv
