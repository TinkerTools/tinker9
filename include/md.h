#ifndef TINKER_MD_H_
#define TINKER_MD_H_

#include "dev_array.h"
#include "energy_buffer.h"
#include "rc_man.h"
#include <string>

TINKER_NAMESPACE_BEGIN
/**
 * \defgroup gvar Global Variables
 *
 * \defgroup md MD Configuration
 * \ingroup gvar
 */

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
/// Padded (by \c MAX_BLOCK_SIZE) number of atoms.
/// \f[ N \le PaddedN \f]
/// \f[ 0 \le PaddedN - N \lt MaxBlockSize \f]
/// \see MAX_BLOCK_SIZE
TINKER_EXTERN int padded_n;

/// \ingroup md
/// \{
/// X, Y, Z coordinates of the current frame or of the entire trajectory frames.
/// \note
/// Arrays are allocated based on \c padded_n.
/// \see padded_n
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
TINKER_EXTERN Energy esum_handle;

/// total potential energy and kinetic energy on host
/// @{
TINKER_EXTERN real epot, eksum, ekin[3][3];

/// total gradients
/// @{
TINKER_EXTERN real *gx, *gy, *gz;
/// @}

/// total virial
TINKER_EXTERN real vir[9];
TINKER_EXTERN Virial vir_handle;

typedef enum {
  thermo_berendsen,
  thermo_bussi,
  thermo_andersen,
  thermo_nose_hoover_chain,
  thermo_null
} thermostat_t;

/// thermostat
TINKER_EXTERN thermostat_t thermostat;

typedef enum {
  baro_berendsen,
  baro_bussi,
  baro_nose_hoover_chain,
  baro_montecarlo,
  baro_null
} barostat_t;

/// barostat
TINKER_EXTERN barostat_t barostat;

namespace calc {
enum {
  xyz = 0x001,  ///< xyz
  vel = 0x002,  ///< velocity
  mass = 0x004, ///< mass
  traj = 0x008, ///< trajectory

  energy = 0x010, ///<  16 energy
  grad = 0x020,   ///<  32 gradient
  virial = 0x040, ///<  64 virial
  analyz = 0x080, ///< 128 analyze

  v0 = energy,                 ///<  16 similar to version 0 tinker energies
  v1 = energy + grad + virial, ///< 112 similar to version 1 tinker energies
  v3 = energy + analyz,        ///< 144 similar to version 3 tinker energies
  v4 = energy + grad,          ///<  48 energy and gradient only
  v5 = grad,                   ///<  32 gradient only
  v6 = grad + virial,          ///<  96 gradient and virial only

  vmask = energy + grad + virial + analyz, ///< version flag mask

  md = 0x100, ///< md calculation
};
}

template <int USE>
void sanity_check() {
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
TINKER_NAMESPACE_END

#endif
