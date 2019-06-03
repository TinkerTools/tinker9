#ifndef TINKER_UTIL_ERROR_CUDART_H_
#define TINKER_UTIL_ERROR_CUDART_H_

#include "error.h"
#include <cuda_runtime.h>

#define m_tinker_get_3rd_arg_(arg1, arg2, arg3, ...) arg3
#define m_tinker_get_4th_arg_(arg1, arg2, arg3, arg4, ...) arg4

#define m_tinker_always_check_cudart_1_(cucall)                                \
  do {                                                                         \
    cudaError_t cures__ = cucall;                                              \
    if (cures__ != cudaSuccess) {                                              \
      print_backtrace();                                                       \
      const char* msg = cudaGetErrorString(cures__);                           \
      std::string m__ =                                                        \
          format(" {} (errno {}) at {}:{}", msg, cures__, __FILE__, __LINE__); \
      throw FatalError(m__);                                                   \
    }                                                                          \
  } while (0)

#define m_tinker_always_check_cudart_2_(cucall, optmsg)                        \
  do {                                                                         \
    cudaError_t cures__ = cucall;                                              \
    if (cures__ != cudaSuccess) {                                              \
      print_backtrace();                                                       \
      const char* msg = cudaGetErrorName(cures__);                             \
      std::string m__ = format(" {} {} (errno {}) at {}:{}", optmsg, msg,      \
                               cures__, __FILE__, __LINE__);                   \
      throw FatalError(m__);                                                   \
    }                                                                          \
  } while (0)

#define m_tinker_always_check_cudart_3_(cucall, res_t, cu_0)                   \
  do {                                                                         \
    res_t cures__ = cucall;                                                    \
    if (cures__ != cu_0) {                                                     \
      print_backtrace();                                                       \
      std::string m__ = format(" errno {} of type {} at {}:{}", cures__,       \
                               TINKER_STR(res_t), __FILE__, __LINE__);         \
      throw FatalError(m__);                                                   \
    }                                                                          \
  } while (0)

#define m_tinker_always_check_cudart_(...)                                     \
  m_tinker_get_4th_arg_(__VA_ARGS__, m_tinker_always_check_cudart_3_,          \
                        m_tinker_always_check_cudart_2_,                       \
                        m_tinker_always_check_cudart_1_)

#if TINKER_DEBUG || defined(TINKER_ALWAYS_CHECK_CUDART)
#  define check_cudart(...)                                                    \
    m_tinker_always_check_cudart_(__VA_ARGS__)(__VA_ARGS__)
#else
#  define check_cudart(cucall, ...) cucall
#endif

#endif
