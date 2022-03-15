#pragma once
#include "mdprec.h"

namespace tinker {
void kinetic_energy(energy_prec& eksum_out, energy_prec (&ekin_out)[3][3], int n,
   const double* mass, const vel_prec* vx, const vel_prec* vy, const vel_prec* vz);
void kinetic_explicit(T_prec& temp_out, energy_prec& eksum_out, energy_prec (&ekin_out)[3][3],
   const vel_prec* vx, const vel_prec* vy, const vel_prec* vz);
void kinetic(T_prec& temp);

/**
 * \ingroup mdpt
 * \brief Berendsen barostat by scaling the coordinates and box dimensions via
 * coupling to an external constant pressure bath. Code for anisotropic pressure
 * coupling was provided by Guido Raos, Dipartimento di Chimica, Politecnico di
 * Milano, Italy.
 *
 * Literature reference:
 *    - <a href="https://doi.org/10.1063/1.448118">
 *    H. J. C. Berendsen, J. P. M. Postma, W. F. van Gunsteren, A. DiNola,
 *    and J. R. Hauk,
 *    "Molecular dynamics with coupling to an external bath",
 *    J. Chem. Phys., 81, 3684-3690 (1984).
 *    </a>
 */

//====================================================================//

void bussi_thermostat(time_prec dt, T_prec temp);
void bussi_thermostat_acc(time_prec dt, T_prec temp);

extern pos_prec *x_pmonte, *y_pmonte, *z_pmonte;
/**
 * \ingroup mdpt
 * \param temp  Current temperature. If the simulation is also isotermal, this
 * argument is not in use. In stead, the target temperature will be used.
 */
void monte_carlo_barostat(energy_prec epot, T_prec temp);
void monte_carlo_barostat_acc(energy_prec epot, T_prec temp);

void berendsen_barostat(time_prec dt);
void berendsen_barostat_acc(time_prec);

//====================================================================//

/**
 * \ingroup mdpt
 * \brief Applies a velocity correction as needed for the Nose-Hoover Chains
 * at the half time step. Not implemented yet.
 *
 * Literature reference:
 *    - <a href="https://doi.org/10.1080/00268979600100761">
 *    G. J. Martyna, M. E. Tuckerman, D. J. Tobias and M. L. Klein,
 *    "Explicit Reversible Integrators for Extended Systems Dynamics",
 *    Molecular Physics, 87, 1117-1157 (1996).
 *    </a>
 */

/**
 * \ingroup mdpt
 * \brief Applies a box size correction as needed for the Monte Carlo barostat
 * at the half time step.
 *
 * Literature reference:
 *    - <a href="https://doi.org/10.1080/00268977200100031">
 *    I. R. McDonald,
 *    "NpT-ensemble Monte Carlo calculations for binary liquid mixtures",
 *    Molecular Physics, 23, 41-58 (1972).
 *    </a>
 */
}
