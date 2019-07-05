#ifndef TINKER_UTIL_ERROR_H_
#define TINKER_UTIL_ERROR_H_

#include "format_print.h"
#include <stdexcept>

TINKER_NAMESPACE_BEGIN
/**
 * @brief
 * Print the calling stack.
 *
 * @param out_stream  A reference to std::ostream object.
 */
void print_backtrace(std::ostream& out_stream = std::cout);

class FatalError : public std::exception {
private:
  std::string msg_;

public:
  FatalError(const char*);
  FatalError(const std::string&);
  FatalError(const FatalError&);
  const char* what() const noexcept override;
};

#define m_tinker_throw(msg)                                                    \
  do {                                                                         \
    print_backtrace();                                                         \
    std::string m_ = format("{} at {}:{}", msg, __FILE__, __LINE__);           \
    throw FatalError(m_);                                                      \
  } while (0)
TINKER_NAMESPACE_END

#endif
