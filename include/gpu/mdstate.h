#ifndef TINKER_GPU_MDSTATE_H_
#define TINKER_GPU_MDSTATE_H_

#include "defines.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
extern int use_data;
extern int n;

extern real* esum;
extern real vir_xx, vir_yx, vir_zx;
extern real vir_xy, vir_yy, vir_zy;
extern real vir_xz, vir_yz, vir_zz;
extern real* mass;

extern real *x, *y, *z;
extern real *vx, *vy, *vz;
extern real *ax, *ay, *az;
extern real *gx, *gy, *gz;

void n_data();
void xyz_data(int op);
void vel_data(int op);
void accel_data(int op);
void mass_data(int op);
void energy_data(int op);

extern box_st* box;
void box_data(int op);
}
TINKER_NAMESPACE_END

#endif
