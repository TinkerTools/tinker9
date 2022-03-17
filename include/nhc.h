#pragma once
#include "md.h"
#include <tinker/detail/bath.hh>

namespace tinker {
/**
 * \brief Maximum length of the NH chain.
 */
constexpr int maxnose = bath::maxnose;

/**
 * \brief Half time-step NHC thermostat for the Velocity-Verlet integrators.
 *
 * Literature reference:
 *    - <a href="https://doi.org/10.1080/00268979600100761">
 *    G. J. Martyna, M. E. Tuckerman, D. J. Tobias and M. L. Klein,
 *    "Explicit Reversible Integrators for Extended Systems Dynamics",
 *    Molecular Physics, 87, 1117-1157 (1996).
 *    </a>
 *
 * \param dt         Full Time-step.
 * \param nnose      Length of the NH chain; no greater than #maxnose.
 * \param vnh        Velocity array of the NH chain.
 * \param qnh        Mass array of the NH chain.
 * \param g0         Degrees of freedom for the kinetic energy.
 * \param f_kin      Function to calculate and returns the pointer to the
 *                   kinetic energy.
 * \param scale_vel  Function to scale the velocities.
 */
void nhc_isot_96(time_prec dt, int nnose, double* vnh, const double* qnh, double g0,
   double* (*f_kin)(), void (*scale_vel)(double));
}
