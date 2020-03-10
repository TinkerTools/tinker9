#pragma once
#include "darray.h"
#include "energy_buffer.h"
#include "rc_man.h"
#include "torque.h"


TINKER_NAMESPACE_BEGIN
/// local frame definitions
enum
{
   pole_none = 0,
   pole_z_only = 1,
   pole_z_then_x = 2,
   pole_bisector = 3,
   pole_z_bisect = 4,
   pole_3_fold = 5
};

/// @brief
/// local axis type and x,y,z-axis defining atoms for each multipole site
struct LocalFrame
{
   int zaxis;  ///< z-axis defining atom, starting from 0
   int xaxis;  ///< x-axis defining atom, starting from 0
   int yaxis;  ///< y-axis defining atom, starting from ONE
   int polaxe; ///< local frame definition
};

// PME: 0, x, y, z, xx, yy, zz, xy, xz, yz
enum
{
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

TINKER_EXTERN real electric, dielec;

TINKER_EXTERN pointer<LocalFrame> zaxis;

/// traceless Cartesian multipoles in the local frame
TINKER_EXTERN pointer<real, mpl_total> pole;
/// traceless Cartesian multipoles in the global frame
TINKER_EXTERN pointer<real, mpl_total> rpole;

/// x, y, and z components of torques on multipole site
/// @{
TINKER_EXTERN pointer<real> trqx, trqy, trqz;
/// @}
/// internal virial Cartesian tensor due to the torques
TINKER_EXTERN virial_buffer vir_trq;

/// direct induced dipole components at each multipole site
TINKER_EXTERN pointer<real, 3> udir;
/// direct induced dipoles in field used for energy terms
TINKER_EXTERN pointer<real, 3> udirp;
/// mutual induced dipole components at each multipole site
TINKER_EXTERN pointer<real, 3> uind;
/// mutual induced dipoles in field used for energy terms
TINKER_EXTERN pointer<real, 3> uinp;

void elec_data(rc_op op);
int use_elec();
bool use_ewald();
void elec_init(int vers);


void chkpole();
void rotpole();
void chkpole_acc();
void rotpole_acc();
TINKER_NAMESPACE_END
