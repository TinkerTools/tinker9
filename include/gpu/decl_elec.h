#ifndef TINKER_GPU_DECL_ELEC_H_
#define TINKER_GPU_DECL_ELEC_H_

#include "decl_real.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
enum {
  elec_coulomb = 1, /// coulomb interaction
  elec_ewald = 2    /// particle mesh ewald summation
};

/// distance at which switching of the potential begins
extern double mpole_switch_cut;
/// distance at which the potential energy goes to zero
extern double mpole_switch_off;

/// local frame definitions
enum {
  pole_none = 0,
  pole_z_only = 1,
  pole_z_then_x = 2,
  pole_bisector = 3,
  pole_z_bisect = 4,
  pole_3_fold = 5
};
typedef struct local_frame_def_st__ {
  int zaxis;  /// z-axis defining atom, starting from 0
  int xaxis;  /// x-axis defining atom, starting from 0
  int yaxis;  /// y-axis defining atom, starting from ONE
  int polaxe; /// local frame definition
} local_frame_t;
extern local_frame_t* zaxis;

// PME: 0, x, y, z, xx, yy, zz, xy, xz, yz
enum {
  mpl_pme_0 = 0,
  mpl_pme_x = 1,
  mpl_pme_y = 2,
  mpl_pme_z = 3,
  mpl_pme_xx = 4,
  mpl_pme_yy = 5,
  mpl_pme_zz = 6,
  mpl_pme_xy = 7,
  mpl_pme_xz = 8,
  mpl_pme_yz = 9,
  mpl_total = 10,
  mpl_pme_yx = mpl_pme_xy,
  mpl_pme_zx = mpl_pme_xz,
  mpl_pme_zy = mpl_pme_yz
};
/// traceless Cartesian multipoles in the local frame
extern real (*pole)[mpl_total];
/// traceless Cartesian multipoles in the global frame
extern real (*rpole)[mpl_total];

/// x, y, and z components of torques on multipole site
extern real *trqx, *trqy, *trqz;
/// internal virial Cartesian tensor due to the torques
extern real* vir_trq;

/// direct induced dipole components at each multipole site
extern real (*udir)[3];
/// direct induced dipoles in field used for energy terms
extern real (*udirp)[3];
/// mutual induced dipole components at each multipole site
extern real (*uind)[3];
/// mutual induced dipoles in field used for energy terms
extern real (*uinp)[3];

/// @return 0 if no multipole electrostatics is involved; otherise, non-zero
int use_elec();

/**
 * @param op  construct or destruct multipole electrostatics
 *
 * has no effect if no multipole electrostatics is involved
 */
void elec_data(int op);

/**
 * initializes the electrostatics calculation
 *
 * Input:
 * @param ver  selects the code path
 *
 * Output:
 * zero torque (if used);
 * zero torque-related virial (if used);
 * call chkpole() and rotpole();
 * if use pme, initialize some pme data structures.
 */
void elec_init(int ver);

/**
 * takes the torque values on a single site defined by a local coordinate frame
 * and converts to Cartesian forces on the original site and sites specifying
 * the local frame.
 *
 * Input:
 * x, y, and z coordinates;
 * x, y, and z torques;
 * local frame definitions;
 * @param ver  selects the code path
 *
 * Output:
 * According to the value of @param ver, it may add add energy gradients, add
 * torque-related virial, or do nothing.
 */
void torque(int ver);
}
TINKER_NAMESPACE_END

#endif
