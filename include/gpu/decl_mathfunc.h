#ifndef TINKER_GPU_DECL_MATHFUNC_H_
#define TINKER_GPU_DECL_MATHFUNC_H_

#include "util/real_mathfunc.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
/**
 * n-dimensional dot product
 * ans = sum (a[i] * b[i]) for i in [1, 2, ..., n]
 *
 * @param gpu_a  device pointer to array a
 * @param gpu_b  device pointer to array b
 * @param cpu_n  number of elements in each array
 * @return       the dot product to the cpu thread
 */
float dotprod(const float* gpu_a, const float* gpu_b, int cpu_n);
double dotprod(const double* gpu_a, const double* gpu_b, int cpu_n);

/**
 * array[i] = scalar * array[i] for i in [1, 2, ..., n]
 *
 * @param gpu_dst  device pointer to the array
 * @param scal     scalar
 * @param nelem    number of elements in the array
 */
void scale_data(float* gpu_dst, float scal, int nelem);
void scale_data(double* gpu_dst, double scal, int nelem);
}
TINKER_NAMESPACE_END

#endif
