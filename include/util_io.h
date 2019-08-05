#ifndef TINKER_UTIL_IO_H_
#define TINKER_UTIL_IO_H_

#include "macro.h"
#include <string>
#include <vector>
TINKER_NAMESPACE_BEGIN
/// @brief
/// plain ascii text stored by lines
class Text : public std::vector<std::string> {
public:
  /// @brief
  /// white space characters
  static constexpr const char* whitespaces = " \t\n\v\f\r";

  /// @return
  /// True if @c ch is one of the white space characters.
  static bool is_ws(char ch);

  /// @return
  /// a c++ string from a c char array
  template <size_t Len>
  static std::string string(const char (&src)[Len]) {
    return std::string(&src[0], &src[0] + Len);
  }

  /// @brief
  /// replace
  /// @{
  static void replace(std::string& src, std::string old, char r);
  static void replace_by_kv(std::string& src, std::string key,
                            std::string value);
  /// @}

  /// @brief
  /// split a string to a vector of string by @c delimiters
  static std::vector<std::string> split(std::string str,
                                        std::string delimiters = whitespaces);

  /// @brief
  /// transform a string to upper or lower case
  /// @{
  static void upcase(std::string&);
  static void lowcase(std::string&);
  /// @}
};
TINKER_NAMESPACE_END

#include <iostream>
#include <sstream>

TINKER_NAMESPACE_BEGIN
/**
 * @brief
 * read ONE value from a/an string/array/io stream;
 * @c arg will not change until the reading successfully exits
 *
 * @return
 * Non-zero if any error happened.
 */
/// @{
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

template <class Arg, size_t Len>
int read_string_1(Arg& arg, const char (&src)[Len]) {
  return read_string_1(arg, src, Len);
};
/// @}

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

//======================================================================
// formatted output

//======================================================================
// fortran string view

TINKER_NAMESPACE_BEGIN
/**
 * @brief
 * References to the non-NULL-terminated Fortran strings.
 * This class provides a few std::string-like methods to handle the Fortran
 * strings in C++ programs.
 */
class FortranStringView {
private:
  char *const b_, *const e_; // [begin, end)
  FortranStringView() = delete;

  /**
   * If dst != src, copy the first_n characters from src to dst;
   * fill dst with blanks if first_n is less than dstlen;
   * dst is NOT terminated by NULL.
   */
  static void copy_with_blank_(char* dst, size_t dstlen, const char* src,
                               size_t first_n);

  /**
   * Compare to string src of len.
   * The shorter string is filled by blanks prior to comparison.
   */
  bool if_eq_(const char* src, size_t len) const;

public:
  //
  // Constructors.
  //
  template <size_t Len>
  FortranStringView(const char (&src)[Len])
      : b_(const_cast<char*>(src))
      , e_(b_ + Len) {}
  FortranStringView(const char* src, size_t len);
  FortranStringView(const char* src);
  FortranStringView(const std::string& src);

  //
  // Assignment operators.
  //
  template <size_t Len>
  FortranStringView& operator=(const char (&src)[Len]) {
    copy_with_blank_(b_, size(), src, Len);
    return *this;
  }
  FortranStringView& operator=(const char* src);
  FortranStringView& operator=(const std::string& src);
  FortranStringView& operator=(const FortranStringView& src);

  //
  // Comparisons.
  //
  template <size_t Len>
  bool operator==(const char (&src)[Len]) const {
    return if_eq_(src, Len);
  }
  bool operator==(const char* src) const;
  bool operator==(const std::string& src) const;
  bool operator==(const FortranStringView& src) const;

  /// @brief  Max number of characters in the string.
  size_t size() const;
  /// @brief  Length of string, ignoring any trailing blanks.
  size_t len_trim() const;
  /// @brief   Trim the trailing blanks.
  /// @return  Trimmed result in std::string.
  std::string trim() const;
  /**
   * @brief
   * Analogous to Fortran str(x:y) syntax.
   *
   * @param begin1  One-based index for the beginning index.
   * @param back1   One-based index for the inclusive ending index.
   * @return        A new FortranStringView object.
   */
  FortranStringView operator()(int begin1, int back1) const;
  /**
   * @brief
   * Analogous to Fortran str(x:) syntax.
   *
   * @param begin1  One-based index for the beginning index.
   * @return        A new FortranStringView object.
   */
  FortranStringView operator()(int begin1) const;
};

typedef FortranStringView fstr_view;
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
void copyin_tinker_arc(const std::string& arcfile, int first1, int last1,
                       int step);
TINKER_NAMESPACE_END

#endif
