#pragma once


#if defined(TINKER_GFORTRAN)
#   define TINKER_MOD(mod, var) __##mod##_MOD_##var
#   define TINKER_RT(rt)        rt##_
#elif defined(TINKER_IFORT)
#   define TINKER_MOD(mod, var) mod##_mp_##var##_
#   define TINKER_RT(rt)        rt##_
#else
#   error We do not know what Fortran compiler you used to compile the Tinker  \
library. You should implement these two macros (TINKER_MOD and TINKER_RT) here \
to mimic its name mangling.
#endif


/**
 * \def TINKER_STR
 * \ingroup macro
 * Convert a predefined macro `s` to a string `"s"`.
 */
#define TINKER_STR(s)   TINKER_STR1_(s)
#define TINKER_STR1_(s) #s


#define TINKER_GET_1ST_ARG(a1, ...)                         a1
#define TINKER_GET_2ND_ARG(a1, a2, ...)                     a2
#define TINKER_GET_3RD_ARG(a1, a2, a3, ...)                 a3
#define TINKER_GET_4TH_ARG(a1, a2, a3, a4, ...)             a4
#define TINKER_GET_5TH_ARG(a1, a2, a3, a4, a5, ...)         a5
#define TINKER_GET_6TH_ARG(a1, a2, a3, a4, a5, a6, ...)     a6
#define TINKER_GET_7TH_ARG(a1, a2, a3, a4, a5, a6, a7, ...) a7


/**
 * \def TINKER_EXTERN_DEFINITION_FILE
 * \ingroup macro
 * Define this macro to true before this header file being included so that
 * the declarations of the "extern" variables will become definitions in the
 * current compilation unit.
 *
 * \def TINKER_EXTERN
 * \ingroup macro
 * In general, macro `TINKER_EXTERN` expands to `extern`, unless macro
 * `TINKER_EXTERN_DEFINITION_FILE` has been predefined to true.
 * This method is useful to declare and define the global variables.
 * \see TINKER_EXTERN_DEFINITION_FILE
 */
#ifndef TINKER_EXTERN_DEFINITION_FILE
#   define TINKER_EXTERN_DEFINITION_FILE 0
#endif
#if TINKER_EXTERN_DEFINITION_FILE
#   define TINKER_EXTERN
#else
#   define TINKER_EXTERN extern
#endif


/**
 * \def TINKER_NAMESPACE
 * \ingroup macro
 *
 * \def TINKER_NAMESPACE_BEGIN
 * \ingroup macro
 *
 * \def TINKER_NAMESPACE_END
 * \ingroup macro
 */
#define TINKER_NAMESPACE       tinker
#define TINKER_NAMESPACE_BEGIN namespace TINKER_NAMESPACE {
#define TINKER_NAMESPACE_END   }


/**
 * \def TINKER_DEBUG
 * \ingroup macro
 * `TINKER_DEBUG` either expands to 0 or 1. It expands to 1 if and only if
 * `DEBUG` is defined and is not defined to 0.
 * `NDEBUG` is the default and it supersedes `DEBUG` should both of them
 * appear. If `DEBUG` is defined to 0, it is equivalent to having `NDEBUG`
 * defined.
 */
#if defined(_DEBUG) && !defined(DEBUG)
#   define DEBUG _DEBUG
#endif
#if defined(_NDEBUG) && !defined(NDEBUG)
#   define NDEBUG _NDEBUG
#endif
#if !defined(NDEBUG) && !defined(DEBUG)
#   define NDEBUG
#   define TINKER_DEBUG 0
#elif defined(NDEBUG)
#   define TINKER_DEBUG 0
#elif defined(DEBUG)
#   define TINKER_DEBUG1_(VAL) VAL##1
#   define TINKER_DEBUG2_(VAL) TINKER_DEBUG1_(VAL)
#   if TINKER_DEBUG2_(DEBUG) == 1
// DEBUG is defined to empty
#      define TINKER_DEBUG 1
#   elif DEBUG != 0
// DEBUG != 0
#      define TINKER_DEBUG 1
#   else
// DEBUG == 0
#      define TINKER_DEBUG 0
#   endif
#else
#   define TINKER_DEBUG 0
#endif


/**
 * \{
 * \def TINKER_DOUBLE_PRECISION
 * \ingroup macro
 * \def TINKER_MIXED_PRECISION
 * \ingroup macro
 * \def TINKER_SINGLE_PRECISION
 * \ingroup macro
 * Only one of the precision macros can be set to 1 and the rest of them will
 * be set to 0.
 *
 * | Macros | real   | mixed  |
 * |--------|--------|--------|
 * | DOUBLE | double | double |
 * | MIXED  | float  | double |
 * | SINGLE | float  | float  |
 * \}
 */
#ifdef TINKER_DOUBLE_PRECISION
#   undef TINKER_DOUBLE_PRECISION
#   define TINKER_DOUBLE_PRECISION 1
#else
#   define TINKER_DOUBLE_PRECISION 0
#endif
#ifdef TINKER_MIXED_PRECISION
#   undef TINKER_MIXED_PRECISION
#   define TINKER_MIXED_PRECISION 1
#else
#   define TINKER_MIXED_PRECISION 0
#endif
#ifdef TINKER_SINGLE_PRECISION
#   undef TINKER_SINGLE_PRECISION
#   define TINKER_SINGLE_PRECISION 1
#else
#   define TINKER_SINGLE_PRECISION 0
#endif
#if (TINKER_DOUBLE_PRECISION + TINKER_MIXED_PRECISION +                        \
     TINKER_SINGLE_PRECISION) != 1
#   error Error in TINKER_?_PRECISION macros detected.
#endif


/**
 * \def TINKER_HOST
 * \ingroup macro
 * Flag for the CPU-only code.
 * \see TINKER_CUDART
 *
 * \def TINKER_CUDART
 * \ingroup macro
 * Flag for the GPU-enabled code.
 * \see TINKER_HOST
 */
#ifndef TINKER_HOST
#   define TINKER_HOST   0
#   define TINKER_CUDART 1
#else
#   undef TINKER_HOST
#   define TINKER_HOST   1
#   define TINKER_CUDART 0
#endif


// C++11
#ifdef __cplusplus
#   if __cplusplus < 201103L
#      error Must enable C++11.
#   endif
#endif


/**
 * \def restrict
 * \ingroup macro
 * Expand to `__restrict__` in the source code.
 */
#define restrict __restrict__


/**
 * \def CONSTEXPR
 * \ingroup macro
 * If possible, use `if CONSTEXPR` to hint at the chances of optimizations.
 */
#if __cplusplus >= 201703L && defined(__cpp_if_constexpr)
#   define CONSTEXPR constexpr
#else
#   define CONSTEXPR
#endif


/**
 * \def MAYBE_UNUSED
 * \ingroup macro
 * Reduce the "unused variable" warnings from the compiler.
 */
#if __has_cpp_attribute(maybe_unused)
#   define MAYBE_UNUSED [[maybe_unused]]
#else
#   define MAYBE_UNUSED __attribute__((unused))
#endif


TINKER_NAMESPACE_BEGIN
#if TINKER_DOUBLE_PRECISION
#   define TINKER_REAL_SIZE  8
#   define TINKER_MIXED_SIZE 8
using real = double;
using mixed = double;
#endif
#if TINKER_MIXED_PRECISION
#   define TINKER_REAL_SIZE  4
#   define TINKER_MIXED_SIZE 8
using real = float;
using mixed = double;
#endif
#if TINKER_SINGLE_PRECISION
#   define TINKER_REAL_SIZE  4
#   define TINKER_MIXED_SIZE 4
using real = float;
using mixed = float;
#endif


using fixed = unsigned long long;
TINKER_NAMESPACE_END


// Suppress Warnings
#if defined(__INTEL_COMPILER)
// #161: unrecognized #pragma
#pragma warning disable 161

#elif defined(__clang__)
#pragma clang diagnostic ignored "-Wextern-c-compat"

#endif
