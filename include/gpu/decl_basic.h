#ifndef TINKER_GPU_DECL_BASIC_H_
#define TINKER_GPU_DECL_BASIC_H_

#include "util/cxx.h"
#include "util/real_mathfunc.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
/// a new type is defined inside the namespace gpu
typedef real_t__ real;

/// resource management operations
typedef enum {
  rc_dealloc = 0x001, /// deallocate device memory
  rc_alloc = 0x002,   /// allocate device memory
  rc_copyin = 0x004,  /// update device data from host memory, or directly
                      /// initialize device data
  rc_copyout = 0x008, /// update host data from device memory
} rc_t;

//======================================================================
void zero_array(real* dst, int nelem);

//======================================================================
// host-device data transfer

// copyin: copy data from host to device
// copyout: data from device to host
void copyin_data(int* dst, const int* src, int nelem);
void copyin_data(real* dst, const double* src, int nelem);
void copyout_data(int* dst, const int* src, int nelem);
void copyout_data(double* dst, const real* src, int nelem);

// copy all src[c][idx0] to dst[c] (c = 0, 1, ..., nelem-1), ie.
// copy all src(idx0+1,f) to dst(f) (f = 1, 2, ..., nelem)
// idx0 = 0, 1, ..., ndim-1
// shape of dst: real dst[nelem], ie. real*(4 or 8) dst(nelem)
// shape of src: double src[nelm][ndim], ie. real*8 src(ndim,nelem)
void copyin_array2(int idx0, int ndim, real* dst, const double* src, int nelem);
void copyout_array2(int idx0, int ndim, double* dst, const real* src,
                    int nelem);

// shape of dst: double dst[nelem][3], ie. real*8 dst(3,nelem)
// shape of src: real src[nelem][3], ie. real*(4 or 8) src(3,nelem)
void copyout_array3(double (*dst)[3], const real (*src)[3], int nelem);
// dst shall be resized inside this function
void copyout_array3(std::vector<std::array<double, 3>>& dst,
                    const real (*src)[3], int nelem);

// transfer data across two device memory address
void copy_array(int* dst, const int* src, int nelem);
void copy_array(real* dst, const real* src, int nelem);
}
TINKER_NAMESPACE_END

#endif
