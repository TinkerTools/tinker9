#pragma once
#include "dev_array.h"
#include "energy_buffer.h"
#include "md_calc.h"
#include "rc_man.h"
#include "time_scale.h"
#include <string>

TINKER_NAMESPACE_BEGIN
/**
 * \ingroup md
 */
extern int rc_flag;

/// \ingroup md
/// Number of the trajectory frames.
TINKER_EXTERN int trajn;

/// \ingroup md
/// Number of atoms.
TINKER_EXTERN int n;
TINKER_EXTERN int padded_n;

/// \ingroup md
/// \{
/// X, Y, Z coordinates of the current frame or of the entire trajectory frames.
TINKER_EXTERN real *trajx, *trajy, *trajz, *x, *y, *z;
/// \}

void goto_frame(int idx0);
void copyin_arc_file(const std::string& arcfile, int first1, int last1,
                     int step);

/// velocities
/// @{
TINKER_EXTERN real *vx, *vy, *vz;
/// @}

/// atomic mass and inversed mass
/// @}
TINKER_EXTERN real *mass, *massinv;
/// @}

/// total potential energy on device
TINKER_EXTERN real esum;
TINKER_EXTERN energy_buffer esum_buf;

/// total potential energy and kinetic energy on host
/// @{
TINKER_EXTERN real epot, eksum, ekin[3][3];

/// total gradients
/// @{
TINKER_EXTERN real *gx, *gy, *gz;
/// @}

/// total virial
TINKER_EXTERN real vir[9];
TINKER_EXTERN virial_buffer vir_buf;

enum class Thermostat
{
   BERENDSEN,
   BUSSI,
   ANDERSEN,
   NOSE_HOOVER_CHAIN,
   NONE
};
constexpr auto BERENDSEN_THERMOSTAT = Thermostat::BERENDSEN;
constexpr auto BUSSI_THERMOSTAT = Thermostat::BUSSI;
constexpr auto ANDERSEN_THERMOSTAT = Thermostat::ANDERSEN;
constexpr auto NOSE_HOOVER_CHAIN_THERMOSTAT = Thermostat::NOSE_HOOVER_CHAIN;
constexpr auto NONE_THERMOSTAT = Thermostat::NONE;

/// thermostat
TINKER_EXTERN Thermostat thermostat;

enum class Barostat
{
   BERENDSEN,
   BUSSI,
   NOSE_HOOVER_CHAIN,
   MONTE_CARLO,
   NONE
};
constexpr auto BERENDSEN_BAROSTAT = Barostat::BERENDSEN;
constexpr auto BUSSI_BAROSTAT = Barostat::BUSSI;
constexpr auto NOSE_HOOVER_CHAIN_BAROSTAT = Barostat::NOSE_HOOVER_CHAIN;
constexpr auto MONTE_CARLO_BAROSTAT = Barostat::MONTE_CARLO;
constexpr auto NONE_BAROSTAT = Barostat::NONE;

/// barostat
TINKER_EXTERN Barostat barostat;

void egv_data(rc_op op);
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
/// @brief
/// zero out global total energy, gradients, and virial on device
/// @{
void zero_egv(int vers);
void zero_egv();
/// @}

void zero_gradient(DMFlag flag, size_t nelem, real* gx, real* gy, real* gz);

/// @brief
/// sum up potential energies and virials on device
void sum_energies(int vers);
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
// mdsave
void mdsave_data(rc_op);

void mdsave_async(int istep, mixed dt);
void mdsave_synchronize();
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
void md_data(rc_op op);

// integrator
void integrate_data(rc_op);

void kinetic(real& temp);
void temper(mixed dt, real& temp);
/**
 * \ingroup md
 * \brief Applies a velocity correction as needed for the Nose-Hoover Chains,
 * and a box size and velocity correction as needed for the Monte Carlo
 * barostat at the half time step.
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
void mdrest(int istep);

void propagate_xyz(mixed dt, bool check_nblist);
/**
 * \brief v += -g/m dt
 */
void propagate_velocity(mixed dt, const real* grx, const real* gry,
                        const real* grz);
/**
 * \brief v += -g/m dt -g2/m dt2
 */
void propagate_velocity2(mixed dt, const real* grx, const real* gry,
                         const real* grz, mixed dt2, const real* grx2,
                         const real* gry2, const real* grz2);
void propagate(int nsteps, mixed dt_ps);

void velocity_verlet(int istep, mixed dt_ps);

void respa_fast_slow(int istep, mixed dt_ps);
const TimeScaleConfig& respa_tsconfig();
constexpr unsigned RESPA_FAST = 1; // 2**0, fast group shall be 0.
constexpr unsigned RESPA_SLOW = 2; // 2**1, slow group shall be 1.
TINKER_EXTERN real *gx1, *gy1, *gz1;
TINKER_EXTERN real *gx2, *gy2, *gz2;
TINKER_NAMESPACE_END


TINKER_NAMESPACE_BEGIN
void kinetic_acc(real&);
void kinetic_cu(real&);

void propagate_xyz_acc(mixed);

void propagate_velocity_acc(mixed, const real*, const real*, const real*);
void propagate_velocity2_acc(mixed, const real*, const real*, const real*,
                             mixed, const real*, const real*, const real*);

void bussi_thermostat(mixed dt, real temp);
void bussi_thermostat_acc(mixed dt, real temp);

void monte_carlo_barostat_update_nb(real epot);
void monte_carlo_barostat_update_nb_acc(real epot);
TINKER_EXTERN real *x_pmonte, *y_pmonte, *z_pmonte;
TINKER_NAMESPACE_END
