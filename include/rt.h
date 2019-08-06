#ifndef TINKER_UTIL_RT_H_
#define TINKER_UTIL_RT_H_

#ifdef TINKER_HOST
#  include "rt_host.h"
#else
#  include "rt_cudart.h"
#endif

#include "gen_unit.h"

TINKER_NAMESPACE_BEGIN
typedef GenericUnit<FFTPlan, 0> FFTPlanUnit;

/// @brief
/// zero-out, deallocate, or allocate bytes on device
/// @{
void zero_bytes(void* ptr, size_t nbytes);
void dealloc_bytes(void* ptr);
void alloc_bytes(void** ptr, size_t nbytes);
template <class T>
void alloc_bytes(T** ptr, size_t nbytes) {
  return alloc_bytes(reinterpret_cast<void**>(ptr), nbytes);
}
/// @}

/// @brief
/// copy @c nbytes from host to device (copyin), from device to host (copyout),
/// or between two device addresses (copy);
/// will block the calling thread
/// @{
void copyin_bytes(void* dst, const void* src, size_t nbytes);
void copyout_bytes(void* dst, const void* src, size_t nbytes);
void copy_bytes(void* dst, const void* src, size_t nbytes);
/// @}

template <>
struct GenericUnitOp<1> {
  struct Dealloc {
    void operator()(void* ptr) { dealloc_bytes(ptr); }
  };

  struct Alloc {
    void operator()(void** ptr, size_t nbytes) { alloc_bytes(ptr, nbytes); }
  };

  struct Copyin {
    void operator()(void* dst, const void* src, size_t nbytes) {
      copyin_bytes(dst, src, nbytes);
    }
  };
};

/// @brief
/// deallocate, allocate, and synchronize the asynchronous stream
/// @{
void dealloc_stream(Stream);
void alloc_stream(Stream*);
void sync_stream(Stream);
/// @}
/// @brief
/// copy between two device addresses without blocking the calling thread
void copy_bytes_async(void* dst, const void* src, size_t nbytes, Stream s);
TINKER_NAMESPACE_END

// fmtlib
#include <ext/fmt/ostream.h>
#include <iostream>
TINKER_NAMESPACE_BEGIN
template <class Out, class Fmt, class... Ts>
void print(Out& out, const Fmt& fmtstr, const Ts&... args) {
  fmt::print(out, fmtstr, args...);
}

template <class Fmt, class... Ts>
std::string format(const Fmt& fmtstr, const Ts&... args) {
  return fmt::format(fmtstr, args...);
}
TINKER_NAMESPACE_END

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
