#ifndef TINKER_GPU_DECL_ELEC_H_
#define TINKER_GPU_DECL_ELEC_H_

#include "decl_real.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
const int elec_coulomb = 1;
const int elec_ewald = 2;

const int pole_none = 0;
const int pole_z_only = 1;
const int pole_z_then_x = 2;
const int pole_bisector = 3;
const int pole_z_bisect = 4;
const int pole_3_fold = 5;
struct local_frame_def_st {
  int zaxis, xaxis, yaxis;
  int polaxe;
};
extern local_frame_def_st* zaxis;

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
  mpl_pme_yx = mpl_pme_xy,
  mpl_pme_zx = mpl_pme_xz,
  mpl_pme_zy = mpl_pme_yz,
};
const int mpl_total = 10;
extern real (*pole)[mpl_total];
extern real (*rpole)[mpl_total];

extern real *trqx, *trqy, *trqz;
extern real* vir_trq;

int use_elec();
void elec_data(int op);

/**
 * @brief
 * Zero torque (if used), torque-related virial (if used), then call chkpole()
 * and rotpole().
 */
void elec_init(int ver);
void torque(int ver);
}
TINKER_NAMESPACE_END

#endif
