#pragma once
#include "tool/io_print.h"
#include <stdexcept>


namespace tinker {
/**
 * \ingroup error
 * Writes out a set of coordinates to a disk file
 * prior to aborting on a serious error.
 */
void prterr();


/**
 * \ingroup error
 * Prints the call stack to a `FILE` pointer.
 */
void print_backtrace(std::FILE* out = stderr);


/**
 * \ingroup error
 * Translates error code to text.
 */
template <class T>
std::string translate_error_code(T error_num);


/**
 * \ingroup error
 * Errors and exceptions that we do not intend to fix or handle.
 */
class FatalError : public std::exception
{
private:
   std::string msg_;


public:
   FatalError(const char* msg)
      : msg_(msg)
   {}
   FatalError(const std::string& msg)
      : msg_(msg)
   {}
   FatalError(const FatalError& e)
      : msg_(e.msg_)
   {}
   const char* what() const noexcept override
   {
      return msg_.c_str();
   }
};


/**
 * \def TINKER_THROW
 * \ingroup error
 * Throws a fatal error message as a `FatalError` exception.
 */
#define TINKER_THROW(msg)                                                      \
   do {                                                                        \
      print_backtrace();                                                       \
      std::string ms_ = msg;                                                   \
      std::string m_ = format("%s at %s:%d", ms_, __FILE__, __LINE__);         \
      throw FatalError(m_);                                                    \
   } while (0)


#define TINKER_ALWAYS_CHECK_RT_1_(call)                                        \
   do {                                                                        \
      auto res_ = call;                                                        \
      if (res_ != 0) {                                                         \
         print_backtrace();                                                    \
         std::string m_;                                                       \
         std::string msg_ = translate_error_code(res_);                        \
         if (msg_ != "")                                                       \
            m_ = format("Errno %d (%s) at %s:%d", res_, msg_, __FILE__,        \
                        __LINE__);                                             \
         else                                                                  \
            m_ = format("Errno %d at %s:%d", res_, __FILE__, __LINE__);        \
         throw FatalError(m_);                                                 \
      }                                                                        \
   } while (0)
#define TINKER_ALWAYS_CHECK_RT_(...)                                           \
   TINKER_GET_2ND_ARG(__VA_ARGS__, TINKER_ALWAYS_CHECK_RT_1_)


/**
 * \def TINKER_ALWAYS_CHECK_RT
 * \ingroup error
 * Defined to `true` in the source code to enable `check_rt` for the release
 * build.
 * \see check_rt
 */
#ifndef TINKER_ALWAYS_CHECK_RT
#   define TINKER_ALWAYS_CHECK_RT 0
#endif


/**
 * \def check_rt
 * \ingroup error
 * It normally does not do extra work other than the function call it captures,
 * unless if either `TINKER_DEBUG` or `TINKER_ALWAYS_CHECK_RT` is `true`,
 * `check_rt()` will check the error code returned by the function call.
 * \see TINKER_DEBUG
 * \see TINKER_ALWAYS_CHECK_RT
 */
#if TINKER_DEBUG || TINKER_ALWAYS_CHECK_RT
#   define check_rt(...) TINKER_ALWAYS_CHECK_RT_(__VA_ARGS__)(__VA_ARGS__)
#else
#   define check_rt(call, ...) call
#endif


/**
 * \def always_check_rt
 * \ingroup error
 * Always checks returned error code and neglects the related macros.
 * \see check_rt
 */
#define always_check_rt(...) TINKER_ALWAYS_CHECK_RT_(__VA_ARGS__)(__VA_ARGS__)
}
