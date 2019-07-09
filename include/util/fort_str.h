#ifndef TINKER_UTIL_FORT_STR_H_
#define TINKER_UTIL_FORT_STR_H_

#include "cxx.h"

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

#endif
