#pragma once
#include "macro.h"

/// \ingroup prec
/// \def TINKER_DETERMINISTIC_FORCE
/// \brief Logical macro for the underlying type of energy gradients.
///    - If `true`, always use fixed-point arithmetic to accumulate energy
///    gradients, regardless of the underlying type of each individual component.
///    - If `false`, use type `real`, which can either be `float` or `double`
///    based on the precision macro.
///
/// \see TINKER_DOUBLE_PRECISION
/// \see TINKER_MIXED_PRECISION
/// \see TINKER_SINGLE_PRECISION
///
/// In general, evaluating energy, forces, etc. twice, we don't expect to get
/// two identical answers, but we may not care as much because the difference
/// is usually negligible.
/// [See [Why is cos(x) !=
/// cos(y)?](https://isocpp.org/wiki/faq/newbie#floating-point-arith2)]
/// Whereas in MD, two simulations with the same initial configurations can
/// easily diverge due to the accumulated difference. If, for whatever reason,
/// you are willing to elongate the process of the inevitable divergence at the
/// cost of slightly slower simulation speed, a more "deterministic" force (using
/// fixed-point arithmetic) can help.
///
/// To accumulate floating-point values via fixed-point arithmetic, we first
/// scale the values by a large integer then sum only the integer part.
/// \code{.cpp}
/// // fixed sum;
/// // FLT is float or double
/// FLT val2 = val * 0x100000000ull;
/// sum += (fixed)((long long)val2);
/// \endcode
/// To get the floating-point sum, we have to cast the integer sum back to
/// floating-point number, then divide it by the same large integer.
/// \code{.cpp}
/// FLT val2 = (FLT)((long long)sum);
/// FLT answer = val2 / 0x100000000ull;
/// \endcode
#ifndef TINKER_DETERMINISTIC_FORCE
#   if TINKER_DOUBLE_PRECISION
#      define TINKER_DETERMINISTIC_FORCE 0
#   else
#      define TINKER_DETERMINISTIC_FORCE 1
#   endif
#endif

namespace tinker {
/// \typedef fixed
/// \ingroup prec
/// \brief 64-bit unsigned integer type for fixed-point arithmetic.
///
/// \typedef real
/// \ingroup prec
/// \brief Floating-point type with lower precision (not higher than #mixed).
/// \see TINKER_MIXED_PRECISION
///
/// \typedef mixed
/// \ingroup prec
/// \brief Floating-point type with higher precision (not lower than #real).
/// \see TINKER_MIXED_PRECISION
using fixed = unsigned long long;
static_assert(sizeof(fixed) == 8, "");

#if TINKER_DOUBLE_PRECISION
#   define TINKER_REAL_SIZE  8
#   define TINKER_MIXED_SIZE 8
using real = double;
using mixed = double;
#endif
#if TINKER_MIXED_PRECISION
#   define TINKER_REAL_SIZE  4
#   define TINKER_MIXED_SIZE 8
using real = float;
using mixed = double;
#endif
#if TINKER_SINGLE_PRECISION
#   define TINKER_REAL_SIZE  4
#   define TINKER_MIXED_SIZE 4
using real = float;
using mixed = float;
#endif

/// \typedef time_prec
/// \ingroup prec
/// \brief Floating-point type for time.
///
/// \typedef T_prec
/// \ingroup prec
/// \brief Floating-point type for temperature.
///
/// \typedef vel_prec
/// \ingroup prec
/// \brief Floating-point type for velocities.
///
/// \typedef pos_prec
/// \ingroup prec
/// \brief Floating-point type for coordinates.
///
/// \typedef e_prec
/// \ingroup prec
/// \brief Recommended floating-point type for the pairwise energy components.
///
/// \typedef v_prec
/// \ingroup prec
/// \brief Recommended floating-point type for the pairwise virial components.
///
/// \typedef g_prec
/// \ingroup prec
/// \brief Recommended floating-point type for the pairwise gradient components.
///
/// \typedef energy_prec
/// \ingroup prec
/// \brief Floating-point type for total energies.
///
/// \typedef virial_prec
/// \ingroup prec
/// \brief Floating-point type for total virials.
///
/// \typedef grad_prec
/// \ingroup prec
/// \brief Floating-point or fixed-point type for the total gradients.
/// \see TINKER_DETERMINISTIC_FORCE
using time_prec = mixed;
using T_prec = mixed;
using vel_prec = mixed;
using pos_prec = mixed;
using e_prec = real;
using v_prec = real;
using g_prec = real;
using energy_prec = mixed;
using virial_prec = mixed;
#if TINKER_DETERMINISTIC_FORCE
using grad_prec = fixed;
#else
using grad_prec = real;
#endif
}
