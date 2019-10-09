#pragma once


#if defined(TINKER_GFORTRAN)
#   define TINKER_MOD(mod, var) __##mod##_MOD_##var
#   define TINKER_RT(rt) rt##_
#elif defined(TINKER_IFORT)
#   define TINKER_MOD(mod, var) mod##_mp_##var##_
#   define TINKER_RT(rt) rt##_
#else
#   error We do not know what Fortran compiler you used to compile the Tinker  \
library. You should implement these two macros (TINKER_MOD and TINKER_RT) here \
to mimic its name mangling. Similarly, you should implement other functions    \
whenever you see an "unknown fortran compiler error".
#endif


/**
 * \def TINKER_EXTERN_DEFINITION_FILE
 * \ingroup macro
 * Define this macro to \c true before this header file being included so that
 * external variable declarations would become variable definitions in the
 * current compilation unit.
 */
/**
 * \def TINKER_EXTERN
 * \ingroup macro
 * In general, macro \c TINKER_EXTERN expands to \c extern, unless macro \c
 * TINKER_EXTERN_DEFINITION_FILE has been predefined to \c true.
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
 * \{
 * \def TINKER_NAMESPACE
 * \ingroup macro
 * \def TINKER_NAMESPACE_BEGIN
 * \ingroup macro
 * \def TINKER_NAMESPACE_END
 * \ingroup macro
 * These namespace macros would possibly improve the source code indentation and
 * make the curly braces less confusing.
 * \}
 */
#define TINKER_NAMESPACE tinker
#define TINKER_NAMESPACE_BEGIN namespace TINKER_NAMESPACE {
#define TINKER_NAMESPACE_END }


/**
 * \def TINKER_DEBUG
 * \ingroup macro
 * \c TINKER_DEBUG either expands to 0 or 1. It expands to 1 if and only if \c
 * DEBUG is defined and is not defined to 0.
 * \c NDEBUG is the default and it supersedes \c DEBUG should both of them
 * appear together. If \c DEBUG is defined to 0, it is equivalent to having \c
 * NDEBUG defined.
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
#   define TINKER_DEBUG_DO_EXPAND_(VAL) VAL##1
#   define TINKER_DEBUG_EXPAND_(VAL) TINKER_DEBUG_DO_EXPAND_(VAL)

#   if TINKER_DEBUG_EXPAND_(DEBUG) == 1
// DEBUG is defined to empty
#      define TINKER_DEBUG 1
#   elif DEBUG != 0
// DEBUG != 0
#      define TINKER_DEBUG 1
#   else
// DEBUG == 0
#      define TINKER_DEBUG 0
#   endif

#   undef TINKER_DEBUG_DO_EXPAND_
#   undef TINKER_DEBUG_EXPAND_
#else
#   define TINKER_DEBUG 0
#endif


/**
 * \{
 * \def TINKER_DOUBLE_PRECISION
 * \ingroup macro
 * \def TINKER_SINGLE_PRECISION
 * \ingroup macro
 * Based on whether \c TINKER_DOUBLE_PRECISION and \c TINKER_SINGLE_PRECISION
 * being pre-defined, these two macros are set to either 0 or 1 as follows
 *
 * | ifdef D | ifdef S |  D  |  S  |
 * |:-------:|:-------:|:---:|:---:|
 * | true    | true    | 0   | 1   |
 * | false   | true    | 0   | 1   |
 * | true    | false   | 1   | 0   |
 * | false   | false   | 0   | 1   |
 * \}
 */
#if defined(TINKER_DOUBLE_PRECISION) && !defined(TINKER_SINGLE_PRECISION)
#   undef TINKER_DOUBLE_PRECISION
#   define TINKER_DOUBLE_PRECISION 1
#   define TINKER_SINGLE_PRECISION 0
#else
#   ifdef TINKER_DOUBLE_PRECISION
#      undef TINKER_DOUBLE_PRECISION
#   endif
#   ifdef TINKER_SINGLE_PRECISION
#      undef TINKER_SINGLE_PRECISION
#   endif
#   define TINKER_DOUBLE_PRECISION 0
#   define TINKER_SINGLE_PRECISION 1
#endif


TINKER_NAMESPACE_BEGIN
/**
 * \typedef real
 * \ingroup util
 * The default floating point type based on the precision macros. Either defined
 * to \c float or \c double,
 *
 * \see TINKER_DOUBLE_PRECISION
 * \see TINKER_SINGLE_PRECISION
 */
#if TINKER_DOUBLE_PRECISION
typedef double real;
#endif
#if TINKER_SINGLE_PRECISION
typedef float real;
#endif


/**
 * \ingroup util
 * The fixed-point type used to accumulate floating-point numbers.
 */
typedef unsigned long long fixed;


/**
 * \ingroup math
 * Constant to convert a number between floating-point and fixed-point
 * representations.
 * \see fixed
 */
constexpr fixed fixed_point = 0x100000000ull;
TINKER_NAMESPACE_END


/**
 * \def TINKER_HOST
 * \ingroup macro
 * If \c true, use standard C++ runtime library on CPU.
 * \see TINKER_CUDART
 *
 * \def TINKER_CUDART
 * \ingroup macro
 * If \c true, use CUDA runtime library on GPU.
 * \see TINKER_HOST
 */
#ifndef TINKER_HOST
#   define TINKER_HOST 0
#   define TINKER_CUDART 1
#else
#   undef TINKER_HOST
#   define TINKER_HOST 1
#   define TINKER_CUDART 0
#endif


// C++11
#ifdef __cplusplus
#   if __cplusplus < 201103L
#      error Must enable C++11.
#   endif
#endif


/**
 * \def if_constexpr
 * \ingroup macro
 * If possible, use `if constexpr` to hint at the chances of optimizations.
 */
#ifdef __cpp_if_constexpr
#   define if_constexpr if constexpr
#else
#   define if_constexpr if
#endif


/**
 * \def RESTRICT
 * \ingroup macro
 * Expand to \c __restrict__ in the source code.
 */
#define RESTRICT __restrict__


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


/**
 * \def CUDA_DEVICE_FUNCTION
 * \ingroup macro
 * Expand to \c __device__ in the CUDA source code.
 */
#ifdef __CUDACC__
#   define CUDA_DEVICE_FUNCTION __device__
#else
#   define CUDA_DEVICE_FUNCTION
#endif


/**
 * \defgroup macro Macros
 *
 * \defgroup gvar Functions and Global Variables
 *
 * \defgroup util Utilities
 * \ingroup gvar
 */
