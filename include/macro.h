#ifndef TINKER_MACRO_H_
#define TINKER_MACRO_H_

// stringification
#define TINKER_STR_IMPL_(s) #s
#define TINKER_STR(s) TINKER_STR_IMPL_(s)

// fortran compiler
#if defined(TINKER_GFORTRAN)
#  define TINKER_MOD(mod, var) __##mod##_MOD_##var
#  define TINKER_RT(rt) rt##_
#elif defined(TINKER_IFORT)
#  define TINKER_MOD(mod, var) mod##_mp_##var##_
#  define TINKER_RT(rt) rt##_
#else
#  error We do not recognize your Fortran compiler. You should implement these \
two macros (TINKER_MOD and TINKER_RT) here to mimic its name mangling. You     \
should also implement other functions whenever you see an                      \
"unknown fortran compiler error".
#endif

// extern
/**
 * @def TINKER_EXTERN
 * @brief
 * In general, macro @c TINKER_EXTERN expands to @c extern, unless macro @c
 * TINKER_EXTERN_DEFINITION_FILE has been predefined, which is useful to declare
 * and define global variables.
 */
#ifdef TINKER_EXTERN_DEFINITION_FILE
#  define TINKER_EXTERN
#else
#  define TINKER_EXTERN extern
#endif

// namespace
#define TINKER_NAMESPACE tinker
#define TINKER_NAMESPACE_BEGIN namespace TINKER_NAMESPACE {
#define TINKER_NAMESPACE_END }

// Debug
#if defined(_DEBUG) && !defined(DEBUG)
#  define DEBUG _DEBUG
#endif
#if defined(_NDEBUG) && !defined(NDEBUG)
#  define NDEBUG _NDEBUG
#endif
/**
 * @def TINKER_DEBUG
 * @brief
 * @c TINKER_DEBUG expands to either 0 or 1. It expands to 1 if and only if @c
 * DEBUG is defined and is not defined to 0.
 * @c NDEBUG is the default and supersedes @c DEBUG.
 * If @c DEBUG is defined to 0, it is equivalent to having @c NDEBUG defined.
 */
#if !defined(NDEBUG) && !defined(DEBUG)
#  define NDEBUG
#  define TINKER_DEBUG 0
#elif defined(NDEBUG)
#  define TINKER_DEBUG 0
#elif defined(DEBUG)
#  define TINKER_DEBUG_DO_EXPAND_(VAL) VAL##1
#  define TINKER_DEBUG_EXPAND_(VAL) TINKER_DEBUG_DO_EXPAND_(VAL)

#  if TINKER_DEBUG_EXPAND_(DEBUG) == 1
// DEBUG is defined to empty
#    define TINKER_DEBUG 1
#  elif DEBUG != 0
// DEBUG != 0
#    define TINKER_DEBUG 1
#  else
// DEBUG == 0
#    define TINKER_DEBUG 0
#  endif

#  undef TINKER_DEBUG_DO_EXPAND_
#  undef TINKER_DEBUG_EXPAND_
#else
#  define TINKER_DEBUG 0
#endif

// precision
#ifndef TINKER_DOUBLE_PRECISION
#  ifndef TINKER_SINGLE_PRECISION
#    define TINKER_SINGLE_PRECISION
#  endif
#endif
TINKER_NAMESPACE_BEGIN
/**
 * @typedef real
 * @brief
 * Floating point type, defined to either @c float or @c double.
 */
#ifdef TINKER_DOUBLE_PRECISION
typedef double real;
#endif
#ifdef TINKER_SINGLE_PRECISION
typedef float real;
#endif
TINKER_NAMESPACE_END

// fixed-point
TINKER_NAMESPACE_BEGIN
typedef unsigned long long FixedPointType;
const FixedPointType fixed_point = 0x100000000ull;
TINKER_NAMESPACE_END

// warp size
TINKER_NAMESPACE_BEGIN
constexpr int WARP_SIZE = 32;
TINKER_NAMESPACE_END

// host vs device
#ifndef TINKER_HOST
#  define TINKER_CUDART
#endif

// compiler features

// if constexpr
#ifdef __cpp_if_constexpr
#  define if_constexpr if constexpr
#else
#  define if_constexpr if
#endif

// maybe unused
#if __has_cpp_attribute(maybe_unused)
#  define MAYBE_UNUSED [[maybe_unused]]
#else
#  define MAYBE_UNUSED __attribute__((unused))
#endif

#endif
