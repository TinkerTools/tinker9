#ifndef TINKER_UTIL_ARRAY_H_
#define TINKER_UTIL_ARRAY_H_

#include "util_cxx.h"

TINKER_NAMESPACE_BEGIN
void zero_array(int* dst, int nelem);
void zero_array(float* dst, int nelem);
void zero_array(double* dst, int nelem);

// copyin: copy data from host to device
// copyout: data from device to host

void copyin_array(int* dst, const int* src, int nelem);
void copyout_array(int* dst, const int* src, int nelem);

void copyin_array(float* dst, const float* src, int nelem);
void copyout_array(float* dst, const float* src, int nelem);

void copyin_array(float* dst, const double* src, int nelem);
void copyout_array(double* dst, const float* src, int nelem);

void copyin_array(double* dst, const double* src, int nelem);
void copyout_array(double* dst, const double* src, int nelem);

// copy all src[c][idx0] to dst[c] (c = 0, 1, ..., nelem-1), i.e.
// copy all src(idx0+1,f) to dst(f) (f = 1, 2, ..., nelem)
// idx0 = 0, 1, ..., ndim-1
// Shape of dst: dst[nelem], i.e. dst(nelem)
// Shape of src: src[nelm][ndim], i.e. src(ndim,nelem)

void copyin_array2(int idx0, int ndim, float* dst, const float* src, int nelem);
void copyout_array2(int idx0, int ndim, float* dst, const float* src,
                    int nelem);

void copyin_array2(int idx0, int ndim, float* dst, const double* src,
                   int nelem);
void copyout_array2(int idx0, int ndim, double* dst, const float* src,
                    int nelem);

void copyin_array2(int idx0, int ndim, double* dst, const double* src,
                   int nelem);
void copyout_array2(int idx0, int ndim, double* dst, const double* src,
                    int nelem);

// Shape of dst: dst[nelem][3], i.e. dst(3,nelem)
// Shape of src: src[nelem][3], i.e. src(3,nelem)

template <class DT, class ST>
void copyout_array3(DT (*dst)[3], const ST (*src)[3], int natom) {
  copyout_array(&dst[0][0], &src[0][0], 3 * natom);
}

// dst shall be resized inside this function
template <class DT, class ST>
void copyout_array3(std::vector<std::array<DT, 3>>& dst, const ST (*src)[3],
                    int natom) {
  dst.resize(natom);
  copyout_array(&dst[0][0], &src[0][0], 3 * natom);
}

// transfer data across two device memory addresses
void copy_array(int* dst, const int* src, int nelem);
void copy_array(float* dst, const float* src, int nelem);
void copy_array(double* dst, const double* src, int nelem);
TINKER_NAMESPACE_END

#endif
