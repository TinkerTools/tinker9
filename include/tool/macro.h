#pragma once

// C++11
#ifdef __cplusplus
#   if __cplusplus < 201103L
#      error Must enable C++11.
#   endif
#endif

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
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
// #20199-D: unrecognized #pragma in device code
#pragma nv_diag_suppress 20199
#endif

//====================================================================//

/// \ingroup cpp_syntax
/// Expands to \c __restrict__, which is a common C++ extension.
#ifdef __cplusplus
#   define restrict __restrict__
#endif

/// \ingroup cpp_syntax
/// `if constexpr` has been added to C++ since C++17.
/// `if CONSTEXPR` expands to `if constexpr` if this feature is supported.
/// Otherwise it expands to `if`.
#if __cplusplus >= 201703L && defined(__cpp_if_constexpr)
#   define CONSTEXPR constexpr
#else
#   define CONSTEXPR
#endif

/// \ingroup cpp_syntax
/// Reduces the "unused variable" warnings from the compiler.
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

/// \ingroup cpp_syntax
/// Converts a predefined macro \c s to a string \c "s".
#define TINKER_STR(s)   TINKER_STR1_(s)
#define TINKER_STR1_(s) #s

#define TINKER_GET_1ST_ARG(a1, ...)                         a1
#define TINKER_GET_2ND_ARG(a1, a2, ...)                     a2
#define TINKER_GET_3RD_ARG(a1, a2, a3, ...)                 a3
#define TINKER_GET_4TH_ARG(a1, a2, a3, a4, ...)             a4
#define TINKER_GET_5TH_ARG(a1, a2, a3, a4, a5, ...)         a5
#define TINKER_GET_6TH_ARG(a1, a2, a3, a4, a5, a6, ...)     a6
#define TINKER_GET_7TH_ARG(a1, a2, a3, a4, a5, a6, a7, ...) a7

/// \def TINKER_EXTERN_DEFINITION_FILE
/// \ingroup cpp_syntax
/// Define this macro to 1 before this header file being included so
/// that the macro \c TINKER_EXTERN will not be expanded to \c extern.
/// \see TINKER_EXTERN
///
/// \def TINKER_EXTERN
/// \ingroup cpp_syntax
/// Expands to \c extern, unless \c TINKER_EXTERN_DEFINITION_FILE has been
/// predefined to 1. This is useful to declare and define the global variables.
/// \see TINKER_EXTERN_DEFINITION_FILE
#ifndef TINKER_EXTERN_DEFINITION_FILE
#   define TINKER_EXTERN_DEFINITION_FILE 0
#endif
#if TINKER_EXTERN_DEFINITION_FILE
#   define TINKER_EXTERN
#else
#   define TINKER_EXTERN extern
#endif

/// \ingroup cpp_syntax
/// Expands to 0 if the macro `NDEBUG` was predefined.
/// Expands to 1 otherwise.
#ifdef NDEBUG
#   define TINKER_DEBUG 0
#else
#   define TINKER_DEBUG 1
#endif

/// \def TINKER9_DIR
/// \ingroup cpp_syntax
/// Path to this source code directory.
#ifndef TINKER9_DIR
#   error TINKER9_DIR is not set.
#else
#   define TINKER9_DIRSTR TINKER_STR(TINKER9_DIR)
#endif

//====================================================================//

/// \def TINKER_CUDART
/// \ingroup platform
/// Macro for the CUDA runtime-enabled GPU code.
/// Defined to 0 when the code is compiled on the CPU platform.
#ifdef TINKER_CUDART
#   undef TINKER_CUDART
#   define TINKER_CUDART 1
#else
#   define TINKER_CUDART 0
#endif

/// \def TINKER_GPULANG_OPENACC
/// \ingroup platform
/// Macro for the OpenACC GPU code.
/// Defined to 0 when
///    - the code is compiled on the CPU platform;
///    - only the CUDA code is in use on the GPU platform.
#ifdef TINKER_GPULANG_OPENACC
#   undef TINKER_GPULANG_OPENACC
#   define TINKER_GPULANG_OPENACC 1
#else
#   define TINKER_GPULANG_OPENACC 0
#endif

/// \def TINKER_GPULANG_CUDA
/// \ingroup platform
/// Macro for the CUDA GPU code.
/// Defined to 0 when
///    - the code is compiled on the CPU platform;
///    - OpenACC code is in use on the GPU platform.
#ifdef TINKER_GPULANG_CUDA
#   undef TINKER_GPULANG_CUDA
#   define TINKER_GPULANG_CUDA 1
#else
#   define TINKER_GPULANG_CUDA 0
#endif

/// \def TINKER_DOUBLE_PRECISION
/// \ingroup prec
/// Macro for the precision mode. Types `real` and `mixed` have different
/// definitions in different modes.
///
/// | Macros | real   | mixed  |
/// |--------|--------|--------|
/// | DOUBLE | double | double |
/// | MIXED  | float  | double |
/// | SINGLE | float  | float  |
///
/// Only one of the precision macros can be set to 1 and the others will be set to 0.
///
/// \def TINKER_MIXED_PRECISION
/// \ingroup prec
/// \copydoc TINKER_DOUBLE_PRECISION
///
/// \def TINKER_SINGLE_PRECISION
/// \ingroup prec
/// \copydoc TINKER_DOUBLE_PRECISION
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
#if (TINKER_DOUBLE_PRECISION + TINKER_MIXED_PRECISION + TINKER_SINGLE_PRECISION) != 1
#   error Detected errors in TINKER_?_PRECISION macros.
#endif
