#ifndef TINKER_GPU_DATA_H_
#define TINKER_GPU_DATA_H_

#include "prec.h"
#include "util/cxx.h"

TINKER_NAMESPACE_BEGIN
const int use_xyz = 0x0001;    /// xyz
const int use_vel = 0x0002;    /// velocity
const int use_accel = 0x0004;  /// acceleration
const int use_mass = 0x0008;   /// mass
const int use_energy = 0x0010; /// energy
const int use_grad = 0x0020;   /// gradient
const int use_virial = 0x0040; /// virial

namespace gpu {
const int op_destroy = 0x00;
const int op_create = 0x01;
void copyin_data(real* dst, const double* src, int nelem);
void copyin_data3(real* d1, real* d2, real* d3, const double* src3, int nelem);

extern int use_data;

extern real *x, *y, *z;
void xyz_data(int op);

extern real *vx, *vy, *vz;
void vel_data(int op);

extern real *ax, *ay, *az;
void accel_data(int op);

extern real* mass;
void mass_data(int op);

extern real esum;

extern real *gx, *gy, *gz;
void grad_data(int op);

extern real vir[3][3];
}

void gpu_data_create();
void gpu_data_destroy();
TINKER_NAMESPACE_END

#endif
