#ifndef TINKER_GPU_DATA_H_
#define TINKER_GPU_DATA_H_

#include "defines.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void copyin_data_1(int* dst, const int* src, int nelem);
void copyin_data_1(real* dst, const double* src, int nelem);
void copyin_data_n(int idx0, int ndim, real* dst, const double* src, int nelem);

void copyout_data_1(int* dst, const int* src, int nelem);
void copyout_data_1(double* dst, const real* src, int nelem);
void copyout_data_n(int idx0, int ndim, double* dst, const real* src,
                    int nelem);

void zero_data(real* dst, int nelem);
}
TINKER_NAMESPACE_END

extern "C" {
void tinker_gpu_data_create();
void tinker_gpu_data_destroy();
}

#endif
