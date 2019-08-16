#ifndef TINKER_MATHFUNC_H_
#define TINKER_MATHFUNC_H_

#include "macro.h"
#include <cmath>
#include <cstdlib>

#define INT_ABS(x) abs(x)

#ifdef TINKER_DOUBLE_PRECISION
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

#ifdef TINKER_SINGLE_PRECISION
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
constexpr real twosix = 1.12246204830937298143;    ///< @f$ \sqrt[6]{2} @f$
constexpr real sqrttwo = 1.41421356237309504880;   ///< @f$ \sqrt{2} @f$
constexpr real sqrtthree = 1.73205080756887729353; ///< @f$ \sqrt{3} @f$

constexpr real elog = M_E;      ///< @f$ exp(1) @f$
constexpr real logten = M_LN10; ///< @f$ ln(10) @f$

constexpr real pi = M_PI;                          ///< @f$ \pi @f$
constexpr real radian = 57.2957795130823208768;    ///< @f$ 180/\pi @f$
constexpr real radinv = 0.01745329251994329576924; ///< @f$ \pi/180 @f$
constexpr real sqrtpi = 1.77245385090551602730;    ///< @f$ \sqrt{\pi} @f$
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
/// @brief
/// find the max or min value of a variadic list
/// @{
template <class T>
T max_of(T a) {
  return a;
}

template <class T, class T2>
T max_of(T a, T2 b) {
  return (a < b) ? b : a;
}

template <class T, class T2, class... Ts>
T max_of(T a, T2 b, Ts... cs) {
  return max_of(max_of(a, b), cs...);
}

template <class T>
T min_of(T a) {
  return a;
}

template <class T, class T2>
T min_of(T a, T2 b) {
  return (a < b) ? a : b;
}

template <class T, class T2, class... Ts>
T min_of(T a, T2 b, Ts... cs) {
  return min_of(min_of(a, b), cs...);
}
/// @}
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
/**
 * @brief
 * n-dimensional dot product
 * @f[
 * dotprod = \sum_i^n a[i] \cdot b[i]
 * @f]
 *
 * @param[in] a
 * device pointer to array a
 *
 * @param[in] b
 * device pointer to array b
 *
 * @param[in] n
 * number of elements in each array
 *
 * @return
 * the dot product to the host thread
 */
/// @{
float dotprod(const float* a, const float* b, int n);
double dotprod(const double* a, const double* b, int n);
/// @}

/**
 * @brief
 * @f[
 * array[i] = scalar \cdot array[i], 1 \leq i \leq n
 * @f]
 *
 * @param[in,out] dst
 * device pointer to the array
 *
 * @param[in] scal
 * scalar
 *
 * @param[in] nelem
 * number of elements in the array
 */
/// @{
void scale_array(float* dst, float scal, int nelem);
void scale_array(double* dst, double scal, int nelem);
/// @}
TINKER_NAMESPACE_END

#endif
