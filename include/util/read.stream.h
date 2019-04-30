#ifndef TINKER_UTIL_READ_STREAM_H_
#define TINKER_UTIL_READ_STREAM_H_

#include "cxx.h"
#include <sstream>

TINKER_NAMESPACE_BEGIN
/**
 * @brief
 * Read ONE value from a string.
 * _arg will not change until the reading successfully exits.
 *
 * @return  Non-zero if any error happened.
 */
template <class __Arg>
int read_string_1(__Arg& _arg, const char* _str, size_t _len) {
  const int succeed = 0;
  const int fail = 1;
  std::istringstream iss(std::string(_str, _len));
  __Arg tmp;
  iss >> tmp;
  if (iss.fail()) {
    return fail;
  } else {
    _arg = tmp;
    return succeed;
  }
}

/**
 * @brief
 * Read one value from a char array.
 * _arg will not change until the reading successfully exits.
 *
 * @return  Non-zero if any error happened.
 */
template <class __Arg, size_t __Len>
int read_string_1(__Arg& _arg, const char (&_src)[__Len]) {
  return read_string_1(_arg, _src, __Len);
};

/**
 * @brief
 * Read one argument from an std::istream object.
 *
 * @param[out] _arg        Argument read from std::istream.
 * @param[in]  _prompt     Prompt string printed on terminal.
 * @param[in]  _auto_fill  Default value to be assigned to _arg.
 * @param[in]  _invalid    Function that returns true if _arg is invalid.
 * @param[in]  _istream    An std::istream object; use std::cin by default.
 */
/*
template <class __Arg, class __Invalid>
void read_stream_1(__Arg& _arg, std::string _prompt, __Arg _auto_fill,
                   __Invalid&& _invalid, std::istream& _istream = std::cin) {
  int input_fail = false;
  std::string line;
  while (_invalid(_arg)) {
    std::cout << _prompt;
    std::getline(_istream, line);
    auto vs = Text::split(line, Text::whitespaces);
    if (vs.size() == 0) {
      _arg = _auto_fill;
    } else {
      input_fail = read_string_1(_arg, line.data(), line.size());
    }
    if (input_fail) {
      // reset failbit
      _istream.clear();
      // expunge the remaining input and invisible '\n'
      // unnecessary if std::getline is used
      // _istream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      goto flag_continue;
    }
    if (_invalid(_arg)) {
      _arg = _auto_fill;
    }
  flag_continue:
    (void)0;
  }
}
*/
TINKER_NAMESPACE_END

#endif
