#ifndef TINKER_ERROR_H_
#define TINKER_ERROR_H_

#include "io_print.h"
#include <stdexcept>

TINKER_NAMESPACE_BEGIN
/// @brief
/// print the call stack
void print_backtrace(std::ostream& out = std::cout);

/// @brief
/// errors and exceptions that we do not intend to fix or handle
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

#define TINKER_THROW(msg)                                                      \
   do {                                                                        \
      print_backtrace();                                                       \
      std::string m_ = format("{} at {}:{}", msg, __FILE__, __LINE__);         \
      throw FatalError(m_);                                                    \
   } while (0)

#define TINKER_GET_3RD_ARG_(arg1, arg2, arg3, ...) arg3

#define TINKER_ALWAYS_CHECK_RT_1_(call)                                        \
   do {                                                                        \
      auto res_ = call;                                                        \
      if (res_ != 0) {                                                         \
         print_backtrace();                                                    \
         std::string m_ =                                                      \
            format(" errno {} at {}:{}", res_, __FILE__, __LINE__);            \
         throw FatalError(m_);                                                 \
      }                                                                        \
   } while (0)

#define TINKER_ALWAYS_CHECK_RT_2_(call, optmsg)                                \
   do {                                                                        \
      auto res_ = call;                                                        \
      if (res_ != 0) {                                                         \
         print_backtrace();                                                    \
         std::string m_ = format(" errno {} ({}) at {}:{}", res_, optmsg,      \
                                 __FILE__, __LINE__);                          \
         throw FatalError(m_);                                                 \
      }                                                                        \
   } while (0)

#define TINKER_ALWAYS_CHECK_RT_(...)                                           \
   TINKER_GET_3RD_ARG_(__VA_ARGS__, TINKER_ALWAYS_CHECK_RT_2_,                 \
                       TINKER_ALWAYS_CHECK_RT_1_)

#if TINKER_DEBUG || defined(TINKER_ALWAYS_CHECK_RT)
#   define check_rt(...) TINKER_ALWAYS_CHECK_RT_(__VA_ARGS__)(__VA_ARGS__)
#else
#   define check_rt(call, ...) call
#endif

#define TINKER_DUMMY_FUNCTION(r)                                               \
   TINKER_THROW(                                                               \
      format("This dummy function, {}(), should not have been called.\n",      \
             TINKER_STR(r)))
TINKER_NAMESPACE_END

#endif
