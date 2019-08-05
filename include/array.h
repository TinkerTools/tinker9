#ifndef TINKER_ARRAY_H_
#define TINKER_ARRAY_H_

#include "macro.h"
#include "rt.h"

TINKER_NAMESPACE_BEGIN
/// @brief
/// zero out the first n elements of an array on device
/// @{
void zero_array(int* dst, int nelem);
void zero_array(float* dst, int nelem);
void zero_array(double* dst, int nelem);
/// @}

// copyin: copy data from host to device
// copyout: data from device to host

/// @brief
/// copy the first n elements bewtween two arrays,
/// either from host to device (copyin),
/// or from device to host (copyout)
/// @{
void copyin_array(int* dst, const int* src, int nelem);
void copyout_array(int* dst, const int* src, int nelem);

void copyin_array(float* dst, const float* src, int nelem);
void copyout_array(float* dst, const float* src, int nelem);

void copyin_array(float* dst, const double* src, int nelem);
void copyout_array(double* dst, const float* src, int nelem);

void copyin_array(double* dst, const double* src, int nelem);
void copyout_array(double* dst, const double* src, int nelem);
/// @}

/// @brief
/// copy the @c idx0-th of every @c ndim elements from @c src to @c dst
///
/// @param[in] idx0
/// ranges from 0 to ndim-1
///
/// @param[in] nelem
/// number of elements copied to @c dst
/// @{
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
/// @}

/// @brief
/// copy the first n elements between two arrays on device
/// @{
void copy_array(int* dst, const int* src, int nelem);
void copy_array(float* dst, const float* src, int nelem);
void copy_array(double* dst, const double* src, int nelem);
/// @}
TINKER_NAMESPACE_END

#endif
