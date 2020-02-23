#pragma once
#include "io_print.h"
#include <stdexcept>


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup error
 * \brief Print the call stack.
 */
void print_backtrace(std::ostream& out = std::cout);


template <class T>
std::string translate_error_message(T error_num);


/**
 * \ingroup error
 * \brief Errors and exceptions that we do not intend to fix or handle.
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
 * Throw a fatal error message as a `FatalError` exception.
 */
#define TINKER_THROW(msg)                                                      \
   do {                                                                        \
      print_backtrace();                                                       \
      std::string m_ = format("{} at {}:{}", msg, __FILE__, __LINE__);         \
      throw FatalError(m_);                                                    \
   } while (0)


#define TINKER_ALWAYS_CHECK_RT_1_(call)                                        \
   do {                                                                        \
      auto res_ = call;                                                        \
      if (res_ != 0) {                                                         \
         print_backtrace();                                                    \
         std::string m_;                                                       \
         std::string msg_ = translate_error_message(res_);                     \
         if (msg_ != "")                                                       \
            m_ = format("Errno {} ({}) at {}:{}", res_, msg_, __FILE__,        \
                        __LINE__);                                             \
         else                                                                  \
            m_ = format("Errno {} at {}:{}", res_, __FILE__, __LINE__);        \
         throw FatalError(m_);                                                 \
      }                                                                        \
   } while (0)
#define TINKER_ALWAYS_CHECK_RT_(...)                                           \
   TINKER_GET_2ND_ARG(__VA_ARGS__, TINKER_ALWAYS_CHECK_RT_1_)

/**
 * \def TINKER_ALWAYS_CHECK_RT
 * \ingroup error
 * Define it to `true` in the source code to enable `check_rt` with the release
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


#define always_check_rt(...) TINKER_ALWAYS_CHECK_RT_(__VA_ARGS__)(__VA_ARGS__)
TINKER_NAMESPACE_END
