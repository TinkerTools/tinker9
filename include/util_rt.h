#ifndef TINKER_UTIL_RT_H_
#define TINKER_UTIL_RT_H_

#ifdef TINKER_HOSTONLY
#  include "util_rt_hostonly.h"
#else
#  include "util_rt_cudart.h"
#endif

#include <stdexcept>
#include "util_io.h"

TINKER_NAMESPACE_BEGIN
/// @brief
/// print the call stack
void print_backtrace(std::ostream& out = std::cout);

/// @brief
/// errors and exceptions that we do not intend to fix
class FatalError : public std::exception {
private:
  std::string msg_;

public:
  FatalError(const char* msg)
      : msg_(msg) {}
  FatalError(const std::string& msg)
      : msg_(msg) {}
  FatalError(const FatalError& e)
      : msg_(e.msg_) {}
  const char* what() const noexcept override { return msg_.c_str(); }
};
TINKER_NAMESPACE_END

#define TINKER_THROW(msg)                                                      \
  do {                                                                         \
    print_backtrace();                                                         \
    std::string m_ = format("{} at {}:{}", msg, __FILE__, __LINE__);           \
    throw FatalError(m_);                                                      \
  } while (0)

#define TINKER_GET_3RD_ARG_(arg1, arg2, arg3, ...) arg3
#define TINKER_GET_4TH_ARG_(arg1, arg2, arg3, arg4, ...) arg4

#define TINKER_ALWAYS_CHECK_CUDART_1_(cucall)                                  \
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

#define TINKER_ALWAYS_CHECK_CUDART_2_(cucall, optmsg)                          \
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

#define TINKER_ALWAYS_CHECK_CUDART_3_(cucall, res_t, cu_0)                     \
  do {                                                                         \
    res_t cures_ = cucall;                                                     \
    if (cures_ != cu_0) {                                                      \
      print_backtrace();                                                       \
      std::string m_ = format(" errno {} of type {} at {}:{}", cures_,         \
                              TINKER_STR(res_t), __FILE__, __LINE__);          \
      throw FatalError(m_);                                                    \
    }                                                                          \
  } while (0)

#define TINKER_ALWAYS_CHECK_CUDART_(...)                                       \
  TINKER_GET_4TH_ARG_(__VA_ARGS__, TINKER_ALWAYS_CHECK_CUDART_3_,              \
                      TINKER_ALWAYS_CHECK_CUDART_2_,                           \
                      TINKER_ALWAYS_CHECK_CUDART_1_)

#if TINKER_DEBUG || defined(TINKER_ALWAYS_CHECK_CUDART)
#  define check_rt(...) TINKER_ALWAYS_CHECK_CUDART_(__VA_ARGS__)(__VA_ARGS__)
#else
#  define check_rt(cucall, ...) cucall
#endif

#endif
