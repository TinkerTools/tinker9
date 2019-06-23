#ifndef TINKER_GPU_DECL_DATAOP_H_
#define TINKER_GPU_DECL_DATAOP_H_

#include "decl_real.h"

/// data operations

TINKER_NAMESPACE_BEGIN
namespace gpu {
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

enum {
  op_dealloc = 0x001, /// deallocate device memory
  op_alloc = 0x002,   /// allocate device memory
  op_copyin = 0x004,  /// update device data from host memory, or directly
                      /// initialize device data
  op_copyout = 0x008, /// update host data from device memory
};
void gpu_data(int op);
}
TINKER_NAMESPACE_END

#endif
