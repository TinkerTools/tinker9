#pragma once
#include "precision.h"

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
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

namespace tinker {
/// \ingroup mdintg
/// \{
/// \var x_pmonte
/// \brief Temporary coordinates created for the Monte Carlo barostat.
/// \var y_pmonte
/// \copydoc x_pmonte
/// \var z_pmonte
/// \copydoc x_pmonte
/// \}
TINKER_EXTERN pos_prec *x_pmonte, *y_pmonte, *z_pmonte;

TINKER_EXTERN double qbar;
TINKER_EXTERN double vbar;
TINKER_EXTERN double vbar_matrix[3][3];
}
