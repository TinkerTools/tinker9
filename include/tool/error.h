#pragma once
#include "tool/ioprint.h"
#include "tool/iotext.h"
#include <stdexcept>

namespace tinker {
/// \ingroup error
/// \brief Writes the current coordinates to a disk file prior to aborting on a serious error.
void printError();

/// \ingroup error
/// \brief Prints the call stack to a `FILE` pointer.
void printBacktrace(std::FILE* out = stderr);

/// \ingroup error
/// \brief Translates the error code to text.
template <class T>
std::string translateErrorCode(T errcode);

/// \ingroup error
/// \brief Errors and exceptions that we do not intend to fix or handle.
class FatalError : public std::exception
{
private:
   std::string m_msg;

public:
   FatalError(const char* msg)
      : m_msg(msg)
   {}
   FatalError(const std::string& msg)
      : m_msg(msg)
   {}
   FatalError(const FatalError& e)
      : m_msg(e.m_msg)
   {}
   const char* what() const noexcept override
   {
      return m_msg.c_str();
   }
};

/// \def TINKER_THROW
/// \ingroup error
/// \brief Throws a fatal error message as a `FatalError` exception.
#define TINKER_THROW(msg)                                                                          \
   do {                                                                                            \
      printBacktrace();                                                                            \
      std::string ms_ = msg;                                                                       \
      std::string m_ = format("%s at %s:%d", ms_, __FILE__, __LINE__);                             \
      throw FatalError(m_);                                                                        \
   } while (0)

#define TINKER_ALWAYS_CHECK_RT_1_(call)                                                            \
   do {                                                                                            \
      auto res_ = call;                                                                            \
      if (res_ != 0) {                                                                             \
         printBacktrace();                                                                         \
         std::string m_;                                                                           \
         std::string msg_ = translateErrorCode(res_);                                              \
         if (msg_ != "")                                                                           \
            m_ = format("Errno %d (%s) at %s:%d", res_, msg_, __FILE__, __LINE__);                 \
         else                                                                                      \
            m_ = format("Errno %d at %s:%d", res_, __FILE__, __LINE__);                            \
         throw FatalError(m_);                                                                     \
      }                                                                                            \
   } while (0)
#define TINKER_ALWAYS_CHECK_RT_(...) TINKER_GET_2ND_ARG(__VA_ARGS__, TINKER_ALWAYS_CHECK_RT_1_)

/// \def TINKER_ALWAYS_CHECK_RT
/// \ingroup error
/// \brief Defined to `true` in the source code to enable `check_rt` for the release build.
/// \see check_rt
#ifndef TINKER_ALWAYS_CHECK_RT
#   define TINKER_ALWAYS_CHECK_RT 0
#endif

/// \def check_rt
/// \ingroup error
/// \brief It normally does not do extra work other than the function call it captures,
/// unless if either `TINKER_DEBUG` or `TINKER_ALWAYS_CHECK_RT` is `true` and this macro
/// will then check the error code returned by the function call.
/// \see TINKER_DEBUG
/// \see TINKER_ALWAYS_CHECK_RT
#if TINKER_DEBUG || TINKER_ALWAYS_CHECK_RT
#   define check_rt(...) TINKER_ALWAYS_CHECK_RT_(__VA_ARGS__)(__VA_ARGS__)
#else
#   define check_rt(call, ...) call
#endif

/// \def always_check_rt
/// \ingroup error
/// \brief Always checks returned error code.
/// \see check_rt
#define always_check_rt(...) TINKER_ALWAYS_CHECK_RT_(__VA_ARGS__)(__VA_ARGS__)
}
