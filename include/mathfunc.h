#ifndef TINKER_MATHFUNC_H_
#define TINKER_MATHFUNC_H_

#include "macro.h"
#include <cmath>
#include <cstdlib>

#define INT_ABS(x) abs(x)

#if TINKER_DOUBLE_PRECISION
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

#if TINKER_SINGLE_PRECISION
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

/// \defgroup math Math
/// \ingroup gvar

TINKER_NAMESPACE_BEGIN
/// \ingroup math
/// \f$ \sqrt[6]{2} \f$
constexpr real twosix = 1.12246204830937298143;
/// \ingroup math
/// \f$ \sqrt{2} \f$
constexpr real sqrttwo = 1.41421356237309504880;
/// \ingroup math
/// \f$ \sqrt{3} \f$
constexpr real sqrtthree = 1.73205080756887729353;

/// \ingroup math
/// \f$ exp(1) \f$
constexpr real elog = M_E;
/// \ingroup math
/// \f$ ln(10) \f$
constexpr real logten = M_LN10;

/// \ingroup math
/// \f$ \pi \f$
constexpr real pi = M_PI;
/// \ingroup math
/// \f$ 180/\pi \f$
constexpr real radian = 57.2957795130823208768;
/// \ingroup math
/// \f$ \pi/180 \f$
constexpr real radinv = 0.01745329251994329576924;
/// \ingroup math
/// \f$ \sqrt{\pi} \f$
constexpr real sqrtpi = 1.77245385090551602730;
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
/// \ingroup math
template <class T>
T max_of(T a) {
  return a;
}

/// \ingroup math
template <class T, class T2>
T max_of(T a, T2 b) {
  return (a < b) ? b : a;
}

/// \ingroup math
/// \return
/// The maximum value from a variadic list of numbers.
template <class T, class T2, class... Ts>
T max_of(T a, T2 b, Ts... cs) {
  return max_of(max_of(a, b), cs...);
}

/// \ingroup math
template <class T>
T min_of(T a) {
  return a;
}

/// \ingroup math
template <class T, class T2>
T min_of(T a, T2 b) {
  return (a < b) ? a : b;
}

/// \ingroup math
/// The minimum value from a variadic list of numbers.
template <class T, class T2, class... Ts>
T min_of(T a, T2 b, Ts... cs) {
  return min_of(min_of(a, b), cs...);
}
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
/// \ingroup math
/// \return
/// True if and only if \c val is a power of 2.
bool is_pow2(size_t val);

/// \ingroup math
/// \return
/// A power of 2 that is less than or equal to \c val.
/// \param val
/// Must be greater than 0.
size_t pow2_le(size_t val);

/// \ingroup math
/// \return
/// A power of 2 that is greater than or equal to \c val.
size_t pow2_ge(size_t val);
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
namespace parallel {
/**
 * \ingroup math
 * \brief Sum all the elements of an array.
 * \return
 * The sum.
 */
int reduce_sum(const int* gpu_a, size_t nelem);
/// \ingroup math
float reduce_sum(const float* gpu_a, size_t nelem);
/// \ingroup math
double reduce_sum(const double* gpu_a, size_t nelem);
/// \ingroup math
unsigned long long reduce_sum(const unsigned long long* gpu_a, size_t nelem);

/**
 * \ingroup math
 * \brief Sum the elements of a 2-dimensional array to an 1-dimensional array.
 * 
 * E.g., a two dimensional array \c v[16][\c m] is used as a virial buffer, in
 * which case, \c nelem = \m, \c neach = 16. The total virial will be written to
 * \c h_ans[\c hn], where \c hn = 9.
 */
void reduce_sum2(int* h_ans, size_t hn, const int* v, size_t nelem,
                 size_t neach);
/// \ingroup math
void reduce_sum2(float* h_ans, size_t hn, const float* v, size_t nelem,
                 size_t neach);
/// \ingroup math
void reduce_sum2(double* h_ans, size_t hn, const double* v, size_t nelem,
                 size_t neach);
/// \ingroup math
void reduce_sum2(unsigned long long* h_ans, size_t hn,
                 const unsigned long long* v, size_t nelem, size_t neach);

/**
 * \ingroup math
 * \brief N-dimensional dot product.
 * \f[ DotProduct = \sum_i^n a_i \cdot b_i \f]
 *
 * \param[in] a
 * Device pointer to array a.
 * \param[in] b
 * Device pointer to array b.
 * \param[in] nelem
 * Number of elements in each array.
 *
 * @return
 * The dot product to the host thread.
 */
float dotprod(const float* a, const float* b, size_t nelem);
/// \ingroup math
double dotprod(const double* a, const double* b, size_t nelem);

/**
 * \ingroup math
 * \brief Multiply all of the elements in an array by a scalar.
 * \f[ a_i = c \cdot a_i, 1 \leq i \leq n \f]
 *
 * \param[in,out] dst
 * Device pointer to the array.
 * \param[in] scal
 * The multiplier, a scalar.
 * \param[in] nelem
 * Number of elements in the array.
 */
void scale_array(float* dst, float scal, size_t nelem);
/// \ingroup math
void scale_array(double* dst, double scal, size_t nelem);
}
TINKER_NAMESPACE_END

#endif
