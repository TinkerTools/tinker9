#pragma once
#include "mdprec.h"


namespace tinker {
void kinetic(T_prec& temp);
void kinetic_acc(T_prec&);
void kinetic_cu(T_prec&);


void temper(time_prec dt, T_prec& temp);


enum class Thermostat
{
   BERENDSEN,
   BUSSI,
   ANDERSEN,
   NOSE_HOOVER_CHAIN,
   LANGEVIN_PISTON,
   NONE
};
constexpr auto BERENDSEN_THERMOSTAT = Thermostat::BERENDSEN;
constexpr auto BUSSI_THERMOSTAT = Thermostat::BUSSI;
constexpr auto ANDERSEN_THERMOSTAT = Thermostat::ANDERSEN;
constexpr auto NOSE_HOOVER_CHAIN_THERMOSTAT = Thermostat::NOSE_HOOVER_CHAIN;
constexpr auto LANGEVIN_PISTON_THERMOSTAT = Thermostat::LANGEVIN_PISTON;
constexpr auto NONE_THERMOSTAT = Thermostat::NONE;
extern Thermostat thermostat;


//====================================================================//


void bussi_thermostat(time_prec dt, T_prec temp);
void bussi_thermostat_acc(time_prec dt, T_prec temp);


enum class Barostat
{
   BERENDSEN,
   BUSSI,
   NOSE_HOOVER_CHAIN,
   LANGEVIN_PISTON,
   MONTE_CARLO,
   NONE
};
constexpr auto BERENDSEN_BAROSTAT = Barostat::BERENDSEN;
constexpr auto BUSSI_BAROSTAT = Barostat::BUSSI;
constexpr auto NOSE_HOOVER_CHAIN_BAROSTAT = Barostat::NOSE_HOOVER_CHAIN;
constexpr auto LANGEVIN_PISTON_BAROSTAT = Barostat::LANGEVIN_PISTON;
constexpr auto MONTE_CARLO_BAROSTAT = Barostat::MONTE_CARLO;
constexpr auto NONE_BAROSTAT = Barostat::NONE;
extern Barostat barostat;


extern pos_prec *x_pmonte, *y_pmonte, *z_pmonte;
extern vel_prec *vx_pmonte, *vy_pmonte, *vz_pmonte;
void monte_carlo_barostat(energy_prec epot);
void monte_carlo_barostat_acc(energy_prec epot);


//====================================================================//


/**
 * \ingroup mdpt
 * \brief Applies a velocity correction as needed for the Nose-Hoover Chains,
 * and a box size correction as needed for the Monte Carlo barostat at the half
 * time step.
 *
 * Literature references:
 *    - G. J. Martyna, M. E. Tuckerman, D. J. Tobias and M. L. Klein,
 *    "Explicit Reversible Integrators for Extended Systems Dynamics",
 *    Molecular Physics, 87, 1117-1157 (1996).
 *    [URL](https://doi.org/10.1080/00268979600100761)
 *    - I. R. McDonald,
 *    "NpT-ensemble Monte Carlo calculations for binary liquid mixtures",
 *    Molecular Physics, 23, 41-58 (1972).
 *    [URL](https://doi.org/10.1080/00268977200100031)
 */
void halftime_correction(bool do_voltrial);
}
