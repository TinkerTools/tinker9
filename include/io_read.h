#ifndef TINKER_IO_READ_H_
#define TINKER_IO_READ_H_

#include "io_text.h"
#include <iostream>
#include <sstream>

TINKER_NAMESPACE_BEGIN
/**
 * @brief
 * read ONE value from a/an string/array and save to @c arg;
 * @c arg will not change until the reading successfully exits
 *
 * @return
 * non-zero if any error happened; otherwise, zero
 */
/// @{
template <class Arg>
int read_string(Arg& arg, const char* str, size_t len)
{
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

template <class Arg, size_t Len>
int read_string(Arg& arg, const char (&src)[Len])
{
   return read_string(arg, src, Len);
};
/// @}

/**
 * @brief
 * @c arg may have been initialized; if its value is invalid, read ONE value
 * from an @c std::istream object; @c arg will not change until the reading
 * successfully exits
 *
 * @param[out] arg
 * values to be read
 *
 * @param[in] prompt
 * message to be sent to @c stdout should input be invalid
 *
 * @param[in] auto_fill
 * default value to be assigned to @c arg
 *
 * @param[in] invalid
 * function that returns @c true if @c arg is invalid or @c false if valid
 *
 * @param[in] istream
 * an @c std::istream object; use std::cin by default
 */
template <class Arg, class Invalid>
void read_stream(Arg& arg, std::string prompt, Arg auto_fill, Invalid&& invalid,
                 std::istream& istream = std::cin)
{
   int input_fail = false;
   std::string line;
   while (invalid(arg)) {
      std::cout << prompt;
      std::getline(istream, line);
      auto vs = Text::split(line);
      if (vs.size() == 0) {
         arg = auto_fill;
      } else {
         input_fail = read_string(arg, line.data(), line.size());
      }
      if (input_fail) {
         // reset failbit
         istream.clear();
         // if using operator<<, expunge the remaining input and invisible '\n';
         // istream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
         // unnecessary if std::getline(...) is used instead
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
