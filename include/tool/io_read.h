#pragma once
#include "tool/io_text.h"
#include <iostream>
#include <sstream>


namespace tinker {
/**
 * \ingroup io
 * Reads ONE value from a string and save to `arg`.
 * `arg` will not change until the reading successfully exits.
 *
 * \return Error code: non-zero if any error occurs; otherwise zero.
 */
/// \{
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
}
/// \}


/**
 * \ingroup io
 * Tests the validity of `arg`. If invalid, read ONE value from an
 * `std::istream` object.
 * `arg` will not change until the reading successfully exits.
 *
 * \param[out] arg       Variable to store the value read from input.
 * \param[in] prompt     Message to be sent to `stdout` should input be invalid.
 * \param[in] auto_fill  Default value to be assigned to `arg`.
 * \param[in] invalid    Function that returns `true` if `arg` is invalid or
 *                       `false` if valid.
 * \param[in] istream    An `std::istream` object; use `std::cin` by default.
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
         istream.clear(); // Reset failbit.
         // Expunge the remaining input and invisible '\n',
         // only if using operator<<;
         // unnecessary if std::getline(...) is used instead.
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
}
