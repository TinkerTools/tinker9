#ifndef TINKER_UTIL_REAL_MATHFUNC_H_
#define TINKER_UTIL_REAL_MATHFUNC_H_

#include "macro.h"
#include <cmath>

TINKER_NAMESPACE_BEGIN

#define INT_ABS(x) abs(x)

#ifdef TINKER_GPU_DOUBLE
typedef double real_t_;

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
typedef float real_t_;

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

constexpr real_t_ twosix = 1.12246204830937298143;    // 2**(1/6)
constexpr real_t_ sqrttwo = 1.41421356237309504880;   // sqrt(2)
constexpr real_t_ sqrtthree = 1.73205080756887729353; // sqrt(3)

constexpr real_t_ elog = M_E;
constexpr real_t_ logten = M_LN10;

constexpr real_t_ pi = M_PI;
constexpr real_t_ radian = 57.2957795130823208768;    // 180/PI
constexpr real_t_ radinv = 0.01745329251994329576924; // PI/180
constexpr real_t_ sqrtpi = 1.77245385090551602730;    // sqrt(PI)

TINKER_NAMESPACE_END

#endif
