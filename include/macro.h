#pragma once


//====================================================================//


#define TINKER_RT(rt) rt##_


//====================================================================//


// C++11
#ifdef __cplusplus
#   if __cplusplus < 201103L
#      error Must enable C++11.
#   endif
#endif


/**
 * \def restrict
 * \ingroup macro
 * Expands to `__restrict__` in the source code.
 */
#define restrict __restrict__


#if defined(__INTEL_COMPILER)
#   define TINKER_ICPC

#elif defined(__PGIC__)
#   define TINKER_PGI

#elif defined(__clang__)
#   define TINKER_CLANG
// xcode clang is different from llvm clang
#   ifdef __apple_build_version__
#      define TINKER_APPLE_CLANG
#   else
#      define TINKER_LLVM_CLANG
#   endif

#elif defined(__GNUC__)
#   define TINKER_GCC
#   if __GNUC__ <= 4 && __GNUC_MINOR__ <= 8
#      warning Your default GNU version is 4.8 where C++11 is incomplete.
#   endif

#endif


// Suppress Warnings
#ifdef TINKER_ICPC
// #161: unrecognized #pragma
#   pragma warning disable 161
#endif
#ifdef TINKER_CLANG
#   pragma clang diagnostic ignored "-Wextern-c-compat"
#endif


/**
 * \def CONSTEXPR
 * \ingroup macro
 * Since C++17, `if constexpr` is added to C++ syntax.
 * `if CONSTEXPR` expands to `if constexpr` if this feature is supported.
 * Otherwise expands to `if`.
 */
#if __cplusplus >= 201703L && defined(__cpp_if_constexpr)
#   define CONSTEXPR constexpr
#else
#   define CONSTEXPR
#endif


/**
 * \def MAYBE_UNUSED
 * \ingroup macro
 * Reduces the "unused variable" warnings from the compiler.
 */
#ifdef __has_cpp_attribute
#   if __has_cpp_attribute(maybe_unused)
#      define MAYBE_UNUSED [[maybe_unused]]
#   else
#      define MAYBE_UNUSED
#   endif
#elif defined(TINKER_ICPC) || defined(TINKER_CLANG) || defined(TINKER_GCC)
#   define MAYBE_UNUSED __attribute__((unused))
#else
#   define MAYBE_UNUSED
#endif


//====================================================================//


/**
 * \def TINKER_STR
 * \ingroup macro
 * Converts a predefined macro `s` to a string `"s"`.
 */
#define TINKER_STR(s)   TINKER_STR1_(s)
#define TINKER_STR1_(s) #s


//====================================================================//


#define TINKER_GET_1ST_ARG(a1, ...)                         a1
#define TINKER_GET_2ND_ARG(a1, a2, ...)                     a2
#define TINKER_GET_3RD_ARG(a1, a2, a3, ...)                 a3
#define TINKER_GET_4TH_ARG(a1, a2, a3, a4, ...)             a4
#define TINKER_GET_5TH_ARG(a1, a2, a3, a4, a5, ...)         a5
#define TINKER_GET_6TH_ARG(a1, a2, a3, a4, a5, a6, ...)     a6
#define TINKER_GET_7TH_ARG(a1, a2, a3, a4, a5, a6, a7, ...) a7


//====================================================================//


/**
 * \def TINKER_EXTERN_DEFINITION_FILE
 * \ingroup macro
 * Define this macro to true before this header files being included so that
 * the declarations of the "extern" variables will become definitions in the
 * current compilation unit.
 *
 * \def TINKER_EXTERN
 * \ingroup macro
 * In general, macro `TINKER_EXTERN` expands to `extern`, unless macro
 * `TINKER_EXTERN_DEFINITION_FILE` has been predefined to true.
 * This is useful to declare and define the global variables.
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


//====================================================================//


/**
 * \def TINKER_DEBUG
 * \ingroup macro
 * Expands to 0 if macro `NDEBUG` was predefined. Expands to 0 otherwise.
 */
#ifdef NDEBUG
#   define TINKER_DEBUG 0
#else
#   define TINKER_DEBUG 1
#endif


//====================================================================//


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


//====================================================================//


/**
 * \def TINKER_DOUBLE_PRECISION
 * \ingroup macro
 * Only one of the precision macros can be set to 1 and the others will
 * be set to 0.
 *
 * | Macros | real   | mixed  |
 * |--------|--------|--------|
 * | DOUBLE | double | double |
 * | MIXED  | float  | double |
 * | SINGLE | float  | float  |
 *
 * \def TINKER_MIXED_PRECISION
 * \ingroup macro
 * \copydoc TINKER_DOUBLE_PRECISION
 *
 * \def TINKER_SINGLE_PRECISION
 * \ingroup macro
 * \copydoc TINKER_DOUBLE_PRECISION
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
#   error Detected errors in TINKER_?_PRECISION macros.
#endif


namespace tinker {
/**
 * \typedef fixed
 * \ingroup prec
 * 64-bit unsigned integer type for fixed-point arithmetic.
 *
 * \typedef real
 * \ingroup prec
 * Floating-point type with lower precision (not higher than #mixed).
 * \see TINKER_MIXED_PRECISION
 *
 * \typedef mixed
 * \ingroup prec
 * Floating-point type with higher precision (not lower than #real).
 * \see TINKER_MIXED_PRECISION
 */
using fixed = unsigned long long;
static_assert(sizeof(fixed) == 8, "");


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
}


//====================================================================//


#ifndef TINKER9_DIR
#   error TINKER9_DIR is not set.
#else
#   define TINKER9_DIRSTR TINKER_STR(TINKER9_DIR)
#endif


//====================================================================//


namespace tinker {
constexpr int MAX_NCHAR = 240;
}
