#ifndef TINKER_GPU_DECL_DATAOP_H_
#define TINKER_GPU_DECL_DATAOP_H_

#include "decl_real.h"

/// data operations

TINKER_NAMESPACE_BEGIN
namespace gpu {
enum {
  op_create = 0x002,
  op_dealloc = 0x001, /// deallocate device memory
  op_alloc = 0x002,   /// allocate device memory
  op_copyin = 0x004   /// update device data from host memory, or directly
                      /// initialize device data
};

void copyin_data(int* dst, const int* src, int nelem);
void copyin_data(real* dst, const double* src, int nelem);
void copyin_data2(int idx0, int ndim, real* dst, const double* src, int nelem);

void copyout_data(int* dst, const int* src, int nelem);
void copyout_data(double* dst, const real* src, int nelem);
void copyout_data2(int idx0, int ndim, double* dst, const real* src, int nelem);
void copyout_data3(double (*dst)[3], const real (*src)[3], int natom);
void copyout_data3(std::vector<std::array<double, 3>>& dst,
                   const real (*src)[3], int natom);

void copy_data(int* dst, const int* src, int nelem);
void copy_data(real* dst, const real* src, int nelem);

void zero_data(real* dst, int nelem);
}
TINKER_NAMESPACE_END

extern "C" {
void tinker_gpu_data_create();
void tinker_gpu_data_destroy();
}

#endif
