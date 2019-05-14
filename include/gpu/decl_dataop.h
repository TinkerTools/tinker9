#ifndef TINKER_GPU_DECL_DATAOP_H_
#define TINKER_GPU_DECL_DATAOP_H_

#include "decl_real.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
const int op_destroy = 0;
const int op_create = 1;

void copyin_data(int* dst, const int* src, int nelem);
void copyin_data(real* dst, const double* src, int nelem);
void copyin_data2(int idx0, int ndim, real* dst, const double* src, int nelem);

void copyout_data(int* dst, const int* src, int nelem);
void copyout_data(double* dst, const real* src, int nelem);
void copyout_data2(int idx0, int ndim, double* dst, const real* src, int nelem);

void zero_data(real* dst, int nelem);

double get_energy(const real* e_gpu);
int get_count(const int* ecount_gpu);
void get_virial(double* v_out, const real* v_gpu);
}
TINKER_NAMESPACE_END

extern "C" {
void tinker_gpu_data_create();
void tinker_gpu_data_destroy();
}

#endif
