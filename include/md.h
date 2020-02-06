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

typedef enum
{
   thermo_berendsen,
   thermo_bussi,
   thermo_andersen,
   thermo_nose_hoover_chain,
   thermo_null
} thermostat_t;

/// thermostat
TINKER_EXTERN thermostat_t thermostat;

typedef enum
{
   baro_berendsen,
   baro_bussi,
   baro_nose_hoover_chain,
   baro_montecarlo,
   baro_null
} barostat_t;

/// barostat
TINKER_EXTERN barostat_t barostat;

template <int USE>
void sanity_check()
{
   constexpr int do_e = USE & calc::energy;
   constexpr int do_a = USE & calc::analyz;
   constexpr int do_g = USE & calc::grad;
   constexpr int do_v = USE & calc::virial;
   // if calc::virial, must calc::grad
   static_assert(do_v ? do_g : true, "");
   // if calc::analyz, must calc::energy
   static_assert(do_a ? do_e : true, "");
}

void egv_data(rc_op op);
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
/// @brief
/// zero out global total energy, gradients, and virial on device
/// @{
void zero_egv(int vers);
void zero_egv();
/// @}

void zero_gradient(int sync, size_t nelem, real* gx, real* gy, real* gz);

/// @brief
/// sum up potential energies and virials on device
void sum_energies(int vers);
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
// mdsave
void mdsave_data(rc_op);

void mdsave_async(int istep, real dt);
void mdsave_synchronize();
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
void md_data(rc_op op);

// integrator
void integrate_data(rc_op);

void kinetic(real& temp);
void temper(real dt, real& temp);
void mdrest(int istep);

void propagate_xyz(real dt, int check_nblist);
/**
 * \brief v += -g/m dt
 */
void propagate_velocity(real dt, const real* grx, const real* gry,
                        const real* grz);
/**
 * \brief v += -g/m dt -g2/m dt2
 */
void propagate_velocity2(real dt, const real* grx, const real* gry,
                         const real* grz, real dt2, const real* grx2,
                         const real* gry2, const real* grz2);
void propagate(int nsteps, real dt_ps);

void velocity_verlet(int istep, real dt_ps);

void respa_fast_slow(int istep, real dt_ps);
const TimeScaleConfig& respa_tsconfig();
constexpr int RESPA_FAST = 1; // 2**0, fast group shall be 0.
constexpr int RESPA_SLOW = 2; // 2**1, slow group shall be 1.
TINKER_EXTERN real *gx1, *gy1, *gz1;
TINKER_EXTERN real *gx2, *gy2, *gz2;

void wait_queue();
TINKER_NAMESPACE_END
