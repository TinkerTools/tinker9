#ifndef TINKER_UTIL_REAL_MATHFUNC_H_
#define TINKER_UTIL_REAL_MATHFUNC_H_

#include "macro.h"
#include <cmath>

TINKER_NAMESPACE_BEGIN

#ifdef TINKER_GPU_DOUBLE
typedef double real_t__;

#  define REAL_SQRT(x) sqrt(x)
#  define REAL_EXP(x) exp(x)
#  define REAL_FLOOR(x) floor(x)
#  define REAL_ABS(x) fabs(x)
#  define REAL_POW(x, idx) pow(x, idx)
#  define REAL_RECIP(x) (1 / ((double)x))
#  define REAL_RSQRT(x) (1 / sqrt((double)x))
#  define REAL_COS(x) cos(x)
#  define REAL_SIN(x) sin(x)
#  define REAL_ERF(x) erf(x)
#  define REAL_ERFC(x) erfc(x)
#  define REAL_MIN(x, y) fmin(x, y)
#  define REAL_MAX(x, y) fmax(x, y)
#endif

#ifdef TINKER_GPU_SINGLE
typedef float real_t__;

#  define REAL_SQRT(x) sqrtf(x)
#  define REAL_EXP(x) expf(x)
#  define REAL_FLOOR(x) floorf(x)
#  define REAL_ABS(x) fabsf(x)
#  define REAL_POW(x, idx) powf(x, idx)
#  define REAL_RECIP(x) (1 / ((float)x))
#  define REAL_RSQRT(x) (1 / sqrtf((float)x))
#  define REAL_COS(x) cosf(x)
#  define REAL_SIN(x) sinf(x)
#  define REAL_ERF(x) erff(x)
#  define REAL_ERFC(x) erfcf(x)
#  define REAL_MIN(x, y) fminf(x, y)
#  define REAL_MAX(x, y) fmaxf(x, y)
#endif

#define INT_ABS(x) abs(x)

#define REAL_SQ(x) ((x) * (x))
#define REAL_CUBE(x) ((x) * (x) * (x))

constexpr real_t__ twosix = 1.12246204830937298143;    // 2**(1/6)
constexpr real_t__ sqrttwo = 1.41421356237309504880;   // sqrt(2)
constexpr real_t__ sqrtthree = 1.73205080756887729353; // sqrt(3)

constexpr real_t__ elog = M_E;
constexpr real_t__ logten = M_LN10;

constexpr real_t__ pi = M_PI;
constexpr real_t__ radian = 57.2957795130823208768; // 180/PI
constexpr real_t__ sqrtpi = 1.77245385090551602730; // sqrt(PI)

TINKER_NAMESPACE_END

#endif
