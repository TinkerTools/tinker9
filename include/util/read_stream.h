#ifndef TINKER_UTIL_READ_STREAM_H_
#define TINKER_UTIL_READ_STREAM_H_

#include "cxx.h"
#include "text.h"
#include <sstream>

TINKER_NAMESPACE_BEGIN
/**
 * @brief
 * Read ONE value from a string.
 * arg_ will not change until the reading successfully exits.
 *
 * @return  Non-zero if any error happened.
 */
template <class Arg_>
int read_string_1(Arg_& arg_, const char* str_, size_t len_) {
  const int succeed = 0;
  const int fail = 1;
  std::istringstream iss(std::string(str_, len_));
  Arg_ tmp;
  iss >> tmp;
  if (iss.fail()) {
    return fail;
  } else {
    arg_ = tmp;
    return succeed;
  }
}

/**
 * @brief
 * Read one value from a char array.
 * arg_ will not change until the reading successfully exits.
 *
 * @return  Non-zero if any error happened.
 */
template <class Arg_, size_t Len_>
int read_string_1(Arg_& arg_, const char (&src_)[Len_]) {
  return read_string_1(arg_, src_, Len_);
};

/**
 * @brief
 * Read one argument from an std::istream object.
 *
 * @param[out] arg_        Argument read from std::istream.
 * @param[in]  prompt_     Prompt string printed on terminal.
 * @param[in]  auto_fill_  Default value to be assigned to arg_.
 * @param[in]  invalid_    Function that returns true if arg_ is invalid.
 * @param[in]  istream_    An std::istream object; use std::cin by default.
 */
template <class Arg_, class Invalid_>
void read_stream_1(Arg_& arg_, std::string prompt_, Arg_ auto_fill_,
                   Invalid_&& invalid_, std::istream& istream_ = std::cin) {
  int input_fail = false;
  std::string line;
  while (invalid_(arg_)) {
    std::cout << prompt_;
    std::getline(istream_, line);
    auto vs = Text::split(line);
    if (vs.size() == 0) {
      arg_ = auto_fill_;
    } else {
      input_fail = read_string_1(arg_, line.data(), line.size());
    }
    if (input_fail) {
      // reset failbit
      istream_.clear();
      // expunge the remaining input and invisible '\n'
      // unnecessary if std::getline is used
      // istream_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      goto flag_continue;
    }
    if (invalid_(arg_)) {
      arg_ = auto_fill_;
    }
  flag_continue:
    (void)0;
  }
}
TINKER_NAMESPACE_END

#endif
