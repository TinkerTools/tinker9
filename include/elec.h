#ifndef TINKER_ELEC_H_
#define TINKER_ELEC_H_

#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
enum {
  elec_coulomb = 1, ///< coulomb interaction
  elec_ewald = 2    ///< particle mesh ewald summation
};

/// local frame definitions
enum {
  pole_none = 0,
  pole_z_only = 1,
  pole_z_then_x = 2,
  pole_bisector = 3,
  pole_z_bisect = 4,
  pole_3_fold = 5
};

struct LocalFrame {
  int zaxis;  ///< z-axis defining atom, starting from 0
  int xaxis;  ///< x-axis defining atom, starting from 0
  int yaxis;  ///< y-axis defining atom, starting from ONE
  int polaxe; ///< local frame definition
};

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

/// distance at which switching of the potential begins
TINKER_EXTERN double mpole_switch_cut;
/// distance at which the potential energy goes to zero
TINKER_EXTERN double mpole_switch_off;

TINKER_EXTERN real electric, dielec;

TINKER_EXTERN LocalFrame* zaxis;

/// traceless Cartesian multipoles in the local frame
TINKER_EXTERN real (*pole)[mpl_total];
/// traceless Cartesian multipoles in the global frame
TINKER_EXTERN real (*rpole)[mpl_total];

/// x, y, and z components of torques on multipole site
/// @{
TINKER_EXTERN real *trqx, *trqy, *trqz;
/// @}
/// internal virial Cartesian tensor due to the torques
TINKER_EXTERN real* vir_trq;

/// direct induced dipole components at each multipole site
TINKER_EXTERN real (*udir)[3];
/// direct induced dipoles in field used for energy terms
TINKER_EXTERN real (*udirp)[3];
/// mutual induced dipole components at each multipole site
TINKER_EXTERN real (*uind)[3];
/// mutual induced dipoles in field used for energy terms
TINKER_EXTERN real (*uinp)[3];

void elec_data(rc_op op);
int use_elec();
void elec_init(int vers);
void torque(int vers);
TINKER_NAMESPACE_END

#endif
