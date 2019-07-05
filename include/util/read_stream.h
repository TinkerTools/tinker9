#ifndef TINKER_UTIL_READ_STREAM_H_
#define TINKER_UTIL_READ_STREAM_H_

#include "cxx.h"
#include "text.h"
#include <sstream>

TINKER_NAMESPACE_BEGIN
/**
 * @brief
 * Read ONE value from a string.
 * arg will not change until the reading successfully exits.
 *
 * @return  Non-zero if any error happened.
 */
template <class Arg>
int read_string_1(Arg& arg, const char* str, size_t len) {
  const int succeed = 0;
  const int fail = 1;
  std::istringstream iss(std::string(str, len));
  Arg tmp;
  iss >> tmp;
  if (iss.fail()) {
    return fail;
  } else {
    arg = tmp;
    return succeed;
  }
}

/**
 * @brief
 * Read one value from a char array.
 * arg will not change until the reading successfully exits.
 *
 * @return  Non-zero if any error happened.
 */
template <class Arg, size_t Len>
int read_string_1(Arg& arg, const char (&src)[Len]) {
  return read_string_1(arg, src, Len);
};

/**
 * @brief
 * Read one argument from an std::istream object.
 *
 * @param[out] arg        Argument read from std::istream.
 * @param[in]  prompt     Prompt string printed on terminal.
 * @param[in]  auto_fill  Default value to be assigned to arg.
 * @param[in]  invalid    Function that returns true if arg is invalid.
 * @param[in]  istream    An std::istream object; use std::cin by default.
 */
template <class Arg, class Invalid>
void read_stream_1(Arg& arg, std::string prompt, Arg auto_fill,
                   Invalid&& invalid, std::istream& istream = std::cin) {
  int input_fail = false;
  std::string line;
  while (invalid(arg)) {
    std::cout << prompt;
    std::getline(istream, line);
    auto vs = Text::split(line);
    if (vs.size() == 0) {
      arg = auto_fill;
    } else {
      input_fail = read_string_1(arg, line.data(), line.size());
    }
    if (input_fail) {
      // reset failbit
      istream.clear();
      // expunge the remaining input and invisible '\n'
      // unnecessary if std::getline is used
      // istream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      goto flag_continue;
    }
    if (invalid(arg)) {
      arg = auto_fill;
    }
  flag_continue:
    (void)0;
  }
}
TINKER_NAMESPACE_END

#endif
