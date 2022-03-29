#pragma once
#include <sstream>
#include <vector>

namespace tinker {
inline namespace v1 {
template <class T>
const T& formatForward(const T& t)
{
   return t;
}

inline const char* formatForward(const std::string& s)
{
   return s.c_str();
}
}

/// \ingroup io_print
/// \brief Formatted output.
/// \param out   A `FILE` pointer.
/// \param f     Format string.
/// \param args  Values to be printed.
template <class F, class... Ts>
void print(std::FILE* out, const F& f, Ts&&... args)
{
   std::fprintf(out, formatForward(f), formatForward(args)...);
}

/// \ingroup io_print
/// \brief Formatted output.
/// \param out  A `FILE` pointer.
/// \param f    Message to be printed.
template <class F>
void print(std::FILE* out, const F& f)
{
   std::fprintf(out, "%s", formatForward(f));
}

/// \ingroup io_print
/// \brief Writes the formatted output to an `std::string` object.
/// \param f     Format string.
/// \param args  Values to be printed out.
/// \return      The formatted `std::string` output.
template <class F, class... Ts>
std::string format(const F& f, Ts&&... args)
{
   const char* fmt = formatForward(f);
   size_t l = std::snprintf(nullptr, 0, fmt, formatForward(args)...);
   size_t sz = l + 1;
   const int len = 128;
   if (sz <= len) {
      char buf[len];
      std::snprintf(buf, sz, fmt, formatForward(args)...);
      return std::string(buf, buf + l);
   } else {
      std::vector<char> ptr(sz);
      char* buf = ptr.data();
      std::snprintf(buf, sz, fmt, formatForward(args)...);
      return std::string(buf, buf + l);
   }
}

/// \ingroup io_print
/// \brief Duplicates then concatenates multiple copies of string.
inline std::string operator*(size_t k, std::string str)
{
   std::ostringstream oss;
   for (size_t i = 0; i < k; ++i) {
      oss << str;
   }
   return oss.str();
}

/// \ingroup io_print
/// \brief Duplicates then concatenate multiple copies of string.
inline std::string operator*(std::string str, size_t k)
{
   return k * str;
}

/// \ingroup io_print
/// \brief Uses `_s` suffix to convert a `const char*` to `std::string`.
inline std::string operator""_s(const char* s, size_t len)
{
   return std::string(s, len);
}
}
