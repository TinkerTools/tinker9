#pragma once
#include "macro.h"
#include <fmt/ostream.h>
#include <iostream>
#include <sstream>
#include <string>


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup io
 * \brief
 * Formatted output using <a href="https://github.com/fmtlib/fmt">fmtlib</a>.
 * \param out    A `FILE*` pointer or an `std::ostream` object.
 * \param fmtstr Format string.
 * \param args   Values to be printed out.
 */
template <class Out, class Fmt, class... Ts>
void print(Out& out, const Fmt& fmtstr, const Ts&... args)
{
   fmt::print(out, fmtstr, args...);
}


/**
 * \ingroup io
 * \brief Write the formatted output to an `std::string` object.
 * \param fmtstr Format string.
 * \param args   Values to be printed out.
 * \return       The formatted `std::string` output.
 */
template <class Fmt, class... Ts>
std::string format(const Fmt& fmtstr, const Ts&... args)
{
   return fmt::format(fmtstr, args...);
}


inline std::string operator*(size_t k, std::string str)
{
   std::ostringstream oss;
   for (size_t i = 0; i < k; ++i) {
      oss << str;
   }
   return oss.str();
}


inline std::string operator*(std::string str, size_t k)
{
   return k * str;
}


inline std::string operator""_s(const char* s, size_t len)
{
   return std::string(s, len);
}
TINKER_NAMESPACE_END
