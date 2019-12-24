#pragma once
#include "dev_array.h"
#include "energy_buffer.h"
#include "rc_man.h"
#include <string>

TINKER_NAMESPACE_BEGIN
/**
 * \ingroup md
 */
TINKER_EXTERN int rc_flag;

/// \ingroup md
/// Number of the trajectory frames.
TINKER_EXTERN int trajn;

/// \ingroup md
/// Number of atoms.
TINKER_EXTERN int n;

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

namespace calc {
/// \ingroup md
/// Use coordinates.
constexpr int xyz = 0x001;
/// \ingroup md
/// Use velocities.
constexpr int vel = 0x002;
/// \ingroup md
/// Use mass.
constexpr int mass = 0x004;
/// \ingroup md
/// Use multi-frame trajectory.
constexpr int traj = 0x008;

/// \ingroup md
/// Evaluate energy.
constexpr int energy = 0x010;
/// \ingroup md
/// Evaluate energy gradient.
constexpr int grad = 0x020;
/// \ingroup md
/// Evaluate virial tensor.
constexpr int virial = 0x040;
/// \ingroup md
/// Evaluate number of interactions.
constexpr int analyz = 0x080;

/// \ingroup md
/// Bits mask to clear energy-irrelevant flags.
constexpr int vmask = energy + grad + virial + analyz;
/// \ingroup md
/// Similar to basic Tinker energy routines.
/// Energy only.
constexpr int v0 = energy;
/// \ingroup md
/// Similar to version 1 Tinker energy routines.
/// Energy, gradient, and virial.
constexpr int v1 = energy + grad + virial;
/// \ingroup md
/// Similar to version 3 Tinker energy routines.
/// Energy and number of interactions.
constexpr int v3 = energy + analyz;
/// \ingroup md
/// Energy and gradient.
constexpr int v4 = energy + grad;
/// \ingroup md
/// Gradient only.
constexpr int v5 = grad;
/// \ingroup md
/// Gradient and virial.
constexpr int v6 = grad + virial;

/// \ingroup md
/// Run MD simulation.
constexpr int md = 0x100;
}

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

void zero_gradient(int nelem, real* gx, real* gy, real* gz);
void zero_gradient_async(int nelem, real* gx, real* gy, real* gz);

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

void propagate_xyz(real dt);
void propagate_velocity(real dt);
void propagate(int nsteps, real dt_ps, void (*itg)(int, real) = nullptr);

void velocity_verlet(int istep, real dt_ps);


void wait_async();
TINKER_NAMESPACE_END
