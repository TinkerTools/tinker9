#ifndef TINKER_UTIL_IO_H_
#define TINKER_UTIL_IO_H_

#include "util_cxx.h"

TINKER_NAMESPACE_BEGIN
class Text : public std::vector<std::string> {
public:
  static constexpr const char* whitespaces = " \t\n\v\f\r";

  /// @return
  /// True if ch is one of the white characters.
  static bool is_ws(char ch);

  template <size_t Len>
  static std::string string(const char (&src)[Len]) {
    return std::string(&src[0], &src[0] + Len);
  }

  // Replace
  static void replace(std::string& s, std::string old, char r);
  static void replace_by_kv(std::string& src, std::string key,
                            std::string value);

  // Split
  static std::vector<std::string> split(std::string str,
                                        std::string delimiters = whitespaces);

  // Case
  static void upcase(std::string&);
  static void lowcase(std::string&);
};

typedef Text text_t;
TINKER_NAMESPACE_END

//======================================================================
// reading from a stream

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

//======================================================================
// formatted output

// fmtlib
#include <ext/fmt/ostream.h>

TINKER_NAMESPACE_BEGIN
template <class Out, class Fmt, class... Ts>
void print(Out& out, const Fmt& fmtstr, const Ts&... args) {
  fmt::print(out, fmtstr, args...);
}

template <class Fmt, class... Ts>
std::string format(const Fmt& fmtstr, const Ts&... args) {
  return fmt::format(fmtstr, args...);
}
TINKER_NAMESPACE_END

//======================================================================
// error

#include <stdexcept>

TINKER_NAMESPACE_BEGIN
/**
 * @brief
 * Print the calling stack.
 *
 * @param out_stream  A reference to std::ostream object.
 */
void print_backtrace(std::ostream& out_stream = std::cout);

class FatalError : public std::exception {
private:
  std::string msg_;

public:
  FatalError(const char* msg) : msg_(msg) {}
  FatalError(const std::string& msg) : msg_(msg) {}
  FatalError(const FatalError& e) : msg_(e.msg_) {}
  const char* what() const noexcept override { return msg_.c_str(); }
};

#define m_tinker_throw(msg)                                                    \
  do {                                                                         \
    print_backtrace();                                                         \
    std::string m_ = format("{} at {}:{}", msg, __FILE__, __LINE__);           \
    throw FatalError(m_);                                                      \
  } while (0)
TINKER_NAMESPACE_END

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
      : b_(const_cast<char*>(src)), e_(b_ + Len) {}
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
namespace gpu {
void copyin_tinker_arc(const std::string& arcfile, int first1, int last1,
                       int step);
}
TINKER_NAMESPACE_END

#endif
