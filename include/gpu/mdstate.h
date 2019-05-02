#ifndef TINKER_GPU_MDSTATE_H_
#define TINKER_GPU_MDSTATE_H_

#include "defines.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void copyin_data_1(int* dst, const int* src, int nelem);
void copyin_data_1(real* dst, const double* src, int nelem);
void copyin_data_n(int idx0, int ndim, real* dst, const double* src, int nelem);

extern int use_data;

extern real* esum;
extern real* vir;
extern real* mass;

extern real *x, *y, *z;
extern real *vx, *vy, *vz;
extern real *ax, *ay, *az;
extern real *gx, *gy, *gz;

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
