#ifndef TINKER_UTIL_MACRO_H_
#define TINKER_UTIL_MACRO_H_

// Stringification
#define TINKER_STR_IMPL_(s) #s
#define TINKER_STR(s) TINKER_STR_IMPL_(s)

// Fortran compiler
#if defined(TINKER_GFORTRAN)
#  define TINKER_MOD(mod, var) __##mod##_MOD_##var
#  define TINKER_RT(rt) rt##_
#elif defined(TINKER_IFORT)
#  define TINKER_MOD(mod, var) mod##_mp_##var##_
#  define TINKER_RT(rt) rt##_
#else
#  error We do not recognize your Fortran compiler. You should implement these \
two macros (TINKER_MOD and TINKER_RT) here to mimic its name mangling.
#endif

// Namespace
#define TINKER_NAMESPACE tinker
#define m_tinker_using_namespace using namespace TINKER_NAMESPACE
#define TINKER_NAMESPACE_BEGIN namespace TINKER_NAMESPACE {
#define TINKER_NAMESPACE_END }

/*
 * Debug
 * If neither DEBUG nor NDEBUG is defined: define NDEBUG; TINKER_DEBUG <- 0.
 * As long as NDEBUG is defined: TINKER_DEBUG <- 0.
 * If DEBUG = 0: define NDEBUG; TINKER_DEBUG <- 0.
 * Otherwise, TINKER_DEBUG <- 1.
 */
#ifdef _DEBUG
#  ifndef DEBUG
#    define DEBUG _DEBUG
#  endif
#endif

#ifdef _NDEBUG
#  ifndef NDEBUG
#    define NDEBUG
#  endif
#endif

#if !defined(NDEBUG) && !defined(DEBUG)
#  define TINKER_DEBUG 0
#elif defined(NDEBUG)
#  define TINKER_DEBUG 0
#elif defined(DEBUG)
#  define TINKER_DEBUG 1
#else
#  error
#endif

#if TINKER_DEBUG == 0 && !defined(NDEBUG)
#  define NDEBUG
#endif

// GPU precision
#ifndef TINKER_GPU_DOUBLE
#  ifndef TINKER_GPU_SINGLE
#    define TINKER_GPU_SINGLE
#  endif
#endif

// features
#ifdef __cpp_if_constexpr
#  define if_constexpr if constexpr
#else
#  define if_constexpr if
#endif

#if __has_cpp_attribute(maybe_unused)
#  define MAYBE_UNUSED [[maybe_unused]]
#else
#  define MAYBE_UNUSED __attribute__((unused))
#endif

#endif
