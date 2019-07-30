#ifndef TINKER_UTIL_MATH_H_
#define TINKER_UTIL_MATH_H_

#include "util_macro.h"
#include <cmath>

#define INT_ABS(x) abs(x)

#ifdef TINKER_GPU_DOUBLE
#  define REAL_SQRT(x) sqrt(x)
#  define REAL_EXP(x) exp(x)
#  define REAL_FLOOR(x) floor(x)
#  define REAL_ABS(x) fabs(x)
#  define REAL_POW(x, idx) pow(x, idx)
#  define REAL_RECIP(x) (1 / ((double)x))
#  define REAL_RSQRT(x) (1 / sqrt((double)x))
#  define REAL_COS(x) cos(x)
#  define REAL_SIN(x) sin(x)
#  define REAL_ACOS(x) acos(x)
#  define REAL_ASIN(x) asin(x)
#  define REAL_ERF(x) erf(x)
#  define REAL_ERFC(x) erfc(x)
#  define REAL_MIN(x, y) fmin(x, y)
#  define REAL_MAX(x, y) fmax(x, y)
#  define REAL_SIGN(x, y) copysign(x, y)
#endif

#ifdef TINKER_GPU_SINGLE
#  define REAL_SQRT(x) sqrtf(x)
#  define REAL_EXP(x) expf(x)
#  define REAL_FLOOR(x) floorf(x)
#  define REAL_ABS(x) fabsf(x)
#  define REAL_POW(x, idx) powf(x, idx)
#  define REAL_RECIP(x) (1 / ((float)x))
#  define REAL_RSQRT(x) (1 / sqrtf((float)x))
#  define REAL_COS(x) cosf(x)
#  define REAL_SIN(x) sinf(x)
#  define REAL_ACOS(x) acosf(x)
#  define REAL_ASIN(x) asinf(x)
#  define REAL_ERF(x) erff(x)
#  define REAL_ERFC(x) erfcf(x)
#  define REAL_MIN(x, y) fminf(x, y)
#  define REAL_MAX(x, y) fmaxf(x, y)
#  define REAL_SIGN(x, y) copysignf(x, y)
#endif

#define REAL_SQ(x) ((x) * (x))
#define REAL_CUBE(x) ((x) * (x) * (x))

TINKER_NAMESPACE_BEGIN
constexpr real twosix = 1.12246204830937298143;    // 2**(1/6)
constexpr real sqrttwo = 1.41421356237309504880;   // sqrt(2)
constexpr real sqrtthree = 1.73205080756887729353; // sqrt(3)

constexpr real elog = M_E;
constexpr real logten = M_LN10;

constexpr real pi = M_PI;
constexpr real radian = 57.2957795130823208768;    // 180/PI
constexpr real radinv = 0.01745329251994329576924; // PI/180
constexpr real sqrtpi = 1.77245385090551602730;    // sqrt(PI)
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
template <class T>
T max_of(T a) {
  return a;
}

template <class T>
T max_of(T a, T b) {
  return (a < b) ? b : a;
}

template <class T, class... Ts>
T max_of(T a, T b, Ts... cs) {
  return max_of(max_of(a, b), cs...);
}

template <class T>
T min_of(T a) {
  return a;
}

template <class T>
T min_of(T a, T b) {
  return (a < b) ? a : b;
}

template <class T, class... Ts>
T min_of(T a, T b, Ts... cs) {
  return min_of(min_of(a, b), cs...);
}
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
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
TINKER_NAMESPACE_END

#endif
