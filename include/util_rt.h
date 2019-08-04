#ifndef TINKER_UTIL_RT_H_
#define TINKER_UTIL_RT_H_

#ifdef TINKER_HOSTONLY
#  include "util_rt_hostonly.h"
#else
#  include "util_rt_cudart.h"
#endif

#include "util_genunit.h"

TINKER_NAMESPACE_BEGIN
typedef GenericUnit<FFTPlan> FFTPlanUnit;

void copyin_bytes(void* dst, const void* src, size_t count);
void copyout_bytes(void* dst, const void* src, size_t count);
void copy_bytes(void* dst, const void* src, size_t count);

void dealloc_stream(Stream);
void alloc_stream(Stream*);
void sync_stream(Stream);
void copy_bytes_async(void* dst, const void* src, size_t count, Stream s);
TINKER_NAMESPACE_END

#include "util_io.h"
#include <stdexcept>

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

#define TINKER_ALWAYS_CHECK_RT_1_(cucall)                                      \
  do {                                                                         \
    auto cures_ = cucall;                                                      \
    if (cures_ != 0) {                                                         \
      print_backtrace();                                                       \
      std::string m_ =                                                         \
          format(" errno {} at {}:{}", cures_, __FILE__, __LINE__);            \
      throw FatalError(m_);                                                    \
    }                                                                          \
  } while (0)

#define TINKER_ALWAYS_CHECK_RT_2_(cucall, optmsg)                              \
  do {                                                                         \
    auto cures_ = cucall;                                                      \
    if (cures_ != 0) {                                                         \
      print_backtrace();                                                       \
      std::string m_ = format(" errno {} ({}) at {}:{}", cures_, optmsg,       \
                              __FILE__, __LINE__);                             \
      throw FatalError(m_);                                                    \
    }                                                                          \
  } while (0)

#define TINKER_ALWAYS_CHECK_RT_(...)                                           \
  TINKER_GET_3RD_ARG_(__VA_ARGS__, TINKER_ALWAYS_CHECK_RT_2_,                  \
                      TINKER_ALWAYS_CHECK_RT_1_)

#if TINKER_DEBUG || defined(TINKER_ALWAYS_CHECK_RT)
#  define check_rt(...) TINKER_ALWAYS_CHECK_RT_(__VA_ARGS__)(__VA_ARGS__)
#else
#  define check_rt(cucall, ...) cucall
#endif

#endif
