#ifndef TINKER_SRC_GPU_RC_H_
#define TINKER_SRC_GPU_RC_H_

#include "rc_cudart.h"
#include "rc_hostonly.h"

//======================================================================

#include "util/error.h"
#include <ext/tinker/tinker_mod.h>
#include <ext/tinker/tinker_rt.h>

#define m_tinker_get_3rd_arg_(arg1, arg2, arg3, ...) arg3
#define m_tinker_get_4th_arg_(arg1, arg2, arg3, arg4, ...) arg4

#define m_tinker_always_check_cudart_1_(cucall)                                \
  do {                                                                         \
    cudaError_t cures_ = cucall;                                               \
    if (cures_ != cudaSuccess) {                                               \
      print_backtrace();                                                       \
      const char* msg = cudaGetErrorString(cures_);                            \
      std::string m_ =                                                         \
          format(" {} (errno {}) at {}:{}", msg, cures_, __FILE__, __LINE__);  \
      throw FatalError(m_);                                                    \
    }                                                                          \
  } while (0)

#define m_tinker_always_check_cudart_2_(cucall, optmsg)                        \
  do {                                                                         \
    cudaError_t cures_ = cucall;                                               \
    if (cures_ != cudaSuccess) {                                               \
      print_backtrace();                                                       \
      const char* msg = cudaGetErrorName(cures_);                              \
      std::string m_ = format(" {} {} (errno {}) at {}:{}", optmsg, msg,       \
                              cures_, __FILE__, __LINE__);                     \
      throw FatalError(m_);                                                    \
    }                                                                          \
  } while (0)

#define m_tinker_always_check_cudart_3_(cucall, res_t, cu_0)                   \
  do {                                                                         \
    res_t cures_ = cucall;                                                     \
    if (cures_ != cu_0) {                                                      \
      print_backtrace();                                                       \
      std::string m_ = format(" errno {} of type {} at {}:{}", cures_,         \
                              TINKER_STR(res_t), __FILE__, __LINE__);          \
      throw FatalError(m_);                                                    \
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
