#ifndef TINKER_ELEC_H_
#define TINKER_ELEC_H_

#include "dev_array.h"
#include "energy_buffer.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
enum class elec_t {
  coulomb, ///< coulomb interaction
  ewald    ///< particle mesh ewald summation
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

/// @brief
/// local axis type and x,y,z-axis defining atoms for each multipole site
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

TINKER_EXTERN device_array::ptr<LocalFrame>::type zaxis;

/// traceless Cartesian multipoles in the local frame
TINKER_EXTERN device_array::ptr<real, mpl_total>::type pole;
/// traceless Cartesian multipoles in the global frame
TINKER_EXTERN device_array::ptr<real, mpl_total>::type rpole;

/// x, y, and z components of torques on multipole site
/// @{
TINKER_EXTERN device_array::ptr<real>::type trqx, trqy, trqz;
/// @}
/// internal virial Cartesian tensor due to the torques
TINKER_EXTERN Virial vir_trq_handle;

/// direct induced dipole components at each multipole site
TINKER_EXTERN device_array::ptr<real, 3>::type udir;
/// direct induced dipoles in field used for energy terms
TINKER_EXTERN device_array::ptr<real, 3>::type udirp;
/// mutual induced dipole components at each multipole site
TINKER_EXTERN device_array::ptr<real, 3>::type uind;
/// mutual induced dipoles in field used for energy terms
TINKER_EXTERN device_array::ptr<real, 3>::type uinp;

void elec_data(rc_op op);
int use_elec();
int use_ewald();
void elec_init(int vers);
void torque(int vers);
TINKER_NAMESPACE_END

#endif
