#pragma once
#include "macro.h"
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace tinker {

//====================================================================//
// FortranStringView

/// \ingroup io
/// Reference to the non-NULL-terminated fortran strings, which
/// provides a few `std::string`-style methods to handle the fortran
/// strings in C++ programs.
class FortranStringView
{
private:
   char* const m_b; // begin in [begin, end)
   char* const m_e; // end in [begin, end)

   FortranStringView() = delete; // no default constructor

   /// If `dst != src`, copy the `first_n` characters from `src` to `dst`;
   /// fill `dst` with blanks if `first_n` is less than `dstlen`;
   /// `dst` is NOT NULL-terminated.
   static void copy_with_blank(char* dst, size_t dstlen, const char* src, size_t first_n);

   /// Compare to a string `src` of `len`.
   /// The shorter string is filled by blanks prior to comparison.
   bool if_eq(const char* src, size_t len) const;

   /// Returns the max number of characters in the string including the trailing whitespaces.
   size_t size() const;

public:
   template <size_t Len>
   FortranStringView(const char (&src)[Len])
      : m_b(const_cast<char*>(src))
      , m_e(m_b + Len)
   {}
   FortranStringView(const char* src, size_t len);
   FortranStringView(const char* src);
   FortranStringView(const std::string& src);

   // assignment
   template <size_t Len>
   FortranStringView& operator=(const char (&src)[Len])
   {
      copy_with_blank(m_b, size(), src, Len);
      return *this;
   }
   FortranStringView& operator=(const char* src);
   FortranStringView& operator=(const std::string& src);
   FortranStringView& operator=(const FortranStringView& src);

   // comparison
   template <size_t Len>
   bool operator==(const char (&src)[Len]) const
   {
      return if_eq(src, Len);
   }
   bool operator==(const char* src) const;
   bool operator==(const std::string& src) const;
   bool operator==(const FortranStringView& src) const;

   /// Returns the length of string, ignoring the trailing whitespaces.
   size_t len_trim() const;

   /// Returns a string without the trailing whitespaces.
   std::string trim() const;

   /// Analogous to Fortran str(x:y) syntax.
   /// \param begin1  One-based index for the beginning index.
   /// \param back1   One-based index for the inclusive ending index.
   /// \return        A new FortranStringView object.
   FortranStringView operator()(int begin1, int back1) const;

   /// Analogous to Fortran str(x:) syntax.
   /// \param begin1  One-based index for the beginning index.
   /// \return        A new FortranStringView object
   FortranStringView operator()(int begin1) const;
};

/// \typedef FstrView
/// \ingroup io
/// An type alias of FortranStringView.
using FstrView = FortranStringView;

//====================================================================//
// print

// Forwards the const reference argument.
template <class T>
const T& formatForward(const T& t)
{
   return t;
}

// Forwards a `const std::string` reference argument to a `const char*`.
inline const char* formatForward(const std::string& s)
{
   return s.c_str();
}

/// \ingroup io
/// Formatted output.
/// \param out   A `FILE` pointer.
/// \param f     Format string.
/// \param args  Values to be printed.
template <class F, class... Ts>
void print(std::FILE* out, const F& f, Ts&&... args)
{
   std::fprintf(out, formatForward(f), formatForward(args)...);
}

/// \ingroup io
/// Formatted output.
/// \param out  A `FILE` pointer.
/// \param f    Message to be printed.
template <class F>
void print(std::FILE* out, const F& f)
{
   std::fprintf(out, "%s", formatForward(f));
}

/// \ingroup io
/// Writes the formatted output to an `std::string` object.
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

/// \ingroup io
/// Duplicates then concatenates multiple copies of string.
inline std::string operator*(size_t k, std::string str)
{
   std::ostringstream oss;
   for (size_t i = 0; i < k; ++i) {
      oss << str;
   }
   return oss.str();
}

/// \ingroup io
/// Duplicates then concatenate multiple copies of string.
inline std::string operator*(std::string str, size_t k)
{
   return k * str;
}

/// \ingroup io
/// Uses `_s` suffix to convert a `const char*` to `std::string`.
inline std::string operator""_s(const char* s, size_t len)
{
   return std::string(s, len);
}

//====================================================================//
// Text

/// \ingroup io
/// Plain ascii text stored by lines.
class Text : public std::vector<std::string>
{
public:
   /// White space characters.
   static constexpr const char* whitespaces = " \t\n\v\f\r";

   /// Returns `true` if `ch` is one of the white space characters.
   static bool is_ws(char ch);

   /// Returns an `std::string` object from a c-style `char` array.
   template <size_t Len>
   static std::string string(const char (&src)[Len])
   {
      return std::string(&src[0], &src[0] + Len);
   }

   /// In string `src`, replaces any characters in `old` to `r`.
   static void replace(std::string& src, std::string old, char r);

   /// In string `src`, replaces any string `key` to string `value`.
   static void replace_by_kv(std::string& src, std::string key, std::string value);

   /// Splits a string `str` by `delimiters` and returns a vector of strings.
   static std::vector<std::string> split(std::string str, std::string delimiters = whitespaces);

   /// Transforms a string to upper case.
   static void upcase(std::string&);

   /// Transforms a string to lower case.
   static void lowcase(std::string&);
};
}

//====================================================================//
// read

#include "tool/io.h"
#include "tool/tinker_suppl.h"

namespace tinker {
/// \ingroup io
/// Reads ONE value from a string and save to `arg`.
/// `arg` will not change until the reading successfully exits.
///
/// \return Error code: non-zero if any error occurs; otherwise zero.
/// \{
template <class Arg>
int ioReadString(Arg& arg, const char* str, size_t len)
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

inline int ioReadString(std::string& arg, const char* str, size_t len)
{
   const int succeed = 0;
   arg = std::string(str, len);
   return succeed;
}

template <class Arg, size_t Len>
int ioReadString(Arg& arg, const char (&src)[Len])
{
   return ioReadString(arg, src, Len);
}
/// \}

/// \ingroup io
/// Tests the validity of `arg`. If invalid, read ONE value from an
/// `std::istream` object.
/// `arg` will not change until the reading successfully exits.
///
/// \param[out] arg       Variable to store the value read from input.
/// \param[in] prompt     Message to be sent to `stdout` should input be invalid.
/// \param[in] auto_fill  Default value to be assigned to `arg`.
/// \param[in] invalid    Function that returns `true` if `arg` is invalid or
///                       `false` if valid.
/// \param[in] istream    An `std::istream` object; use `std::cin` by default.
template <class Arg, class Invalid>
void ioReadStream(
   Arg& arg, std::string prompt, Arg auto_fill, Invalid&& invalid, std::istream& istream = std::cin)
{
   bool is_cin = (&istream == &std::cin);
   int input_fail = invalid(arg); // True if arg is out of range.
   std::string line;
   while (input_fail) {
      std::printf("%s", prompt.c_str());
      std::fflush(stdout);
      if (is_cin) {
         line = tinker_f_read_stdin_line();
      } else {
         istream.clear(); // Reset failbit.
         // Expunge the remaining input and invisible '\n',
         // only if using operator<<;
         // unnecessary if std::getline(...) is used instead.
         // istream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
         std::getline(istream, line);
      }
      auto vs = Text::split(line);
      if (vs.size() == 0) {
         arg = auto_fill;
         input_fail = false;
      } else {
         // True if failed to parse input.
         input_fail = ioReadString(arg, line.data(), line.size());
      }
      if (!input_fail)
         input_fail = invalid(arg);
   }
}

/// \ingroup io
/// Rewind a stream object.
template <class T>
void ioRewindStream(T& stream)
{
   stream.clear();
   stream.seekg(0);
}
}
