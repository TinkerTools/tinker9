#pragma once
#include "macro.h"
#include <fmt/ostream.h>
#include <iostream>
#include <sstream>
#include <string>


TINKER_NAMESPACE_BEGIN
template <class T>
const T& fmtcvt(const T& t)
{
   return t;
}


inline const char* fmtcvt(const std::string& s)
{
   return s.c_str();
}


/**
 * \ingroup io
 * \brief
 * Formatted output using <a href="https://github.com/fmtlib/fmt">fmtlib</a>.
 * \param out    A `FILE*` pointer or an `std::ostream` object.
 * \param fmtstr Format string.
 * \param args   Values to be printed out.
 */
template <class Out, class Fmt, class... Ts>
void PRINT(Out& out, const Fmt& fmtstr, const Ts&... args)
{
   fmt::print(out, fmtstr, args...);
}
template <class F, class... Ts>
void print(std::FILE* out, const F& f, Ts&&... args)
{
   std::fprintf(out, fmtcvt(f), fmtcvt(args)...);
}
template <class F>
void print(std::FILE* out, const F& f)
{
   std::fprintf(out, "%s", fmtcvt(f));
}


/**
 * \ingroup io
 * \brief Write the formatted output to an `std::string` object.
 * \param fmtstr Format string.
 * \param args   Values to be printed out.
 * \return       The formatted `std::string` output.
 */
template <class Fmt, class... Ts>
std::string FORMAT(const Fmt& fmtstr, const Ts&... args)
{
   return fmt::format(fmtstr, args...);
}
template <class F, class... Ts>
std::string format(const F& f, Ts&&... args)
{
   const char* fmt = fmtcvt(f);
   size_t l = std::snprintf(nullptr, 0, fmt, fmtcvt(args)...);
   size_t sz = l + 1;
   const int len = 64;
   if (sz <= len) {
      char buf[len];
      std::snprintf(buf, sz, fmt, fmtcvt(args)...);
      return std::string(buf, buf + l);
   } else {
      std::vector<char> ptr(sz);
      char* buf = ptr.data();
      std::snprintf(buf, sz, fmt, fmtcvt(args)...);
      return std::string(buf, buf + l);
   }
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
