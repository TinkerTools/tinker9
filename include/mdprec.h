#pragma once
#include "macro.h"


/**
 * \ingroup mdprec
 * \page mdprec
 *
 * | Properties                  | Types       | Underlying Types |
 * |----------------------------:|------------:|-----------------:|
 * | time                        | time_prec   | mixed            |
 * | temperature                 | T_prec      | mixed            |
 * | mass                        | mass_prec   | real             |
 * | velocity                    | vel_prec    | mixed            |
 * | position (for integrators)  | pos_prec    | mixed            |
 * | energy (for integrators)    | energy_prec | mixed            |
 * | virial (for integrators)    | virial_prec | mixed            |
 * | gradient (deterministic)    | grad_prec   | fixed            |
 * | gradient (stochastic error) | grad_prec   | real             |
 * | energy (individual terms)   | e_prec      | real             |
 * | virial (individual terms)   | v_prec      | real             |
 * | gradient (individual terms) | undefined   | real             |
 * | position (for energies)     | undefined   | real             |
 */


/**
 * \ingroup mdprec
 * \def TINKER_DETERMINISTIC_FORCE
 * Logical macro for the underlying type of energy gradients.
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
 * [See [Why is cos(x) != cos(y)?](https://isocpp.org/wiki/faq/newbie)]
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
using mass_prec = real;
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
