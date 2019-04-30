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
   * If _dst != _src, copy the _first_n characters from _src to _dst;
   * fill _dst with blanks if _first_n is less than _dstlen;
   * _dst is NOT terminated by NULL.
   */
  static void copy_with_blank_(char* _dst, size_t _dstlen, const char* _src,
                               size_t _first_n);

  /**
   * Compare to string _src of _len.
   * The shorter string is filled by blanks prior to comparison.
   */
  bool if_eq_(const char* _src, size_t _len) const;

public:
  //
  // Constructors.
  //
  template <size_t __Len>
  FortranStringView(const char (&_src)[__Len])
      : b_(const_cast<char*>(_src)), e_(b_ + __Len) {}
  FortranStringView(const char* _src, size_t _len);
  FortranStringView(const char* _src);
  FortranStringView(const std::string& _src);

  //
  // Assignment operators.
  //
  template <size_t __Len>
  FortranStringView& operator=(const char (&_src)[__Len]) {
    copy_with_blank_(b_, size(), _src, __Len);
    return *this;
  }
  FortranStringView& operator=(const char* _src);
  FortranStringView& operator=(const std::string& _src);
  FortranStringView& operator=(const FortranStringView& _src);

  //
  // Comparisons.
  //
  template <size_t __Len>
  bool operator==(const char (&_src)[__Len]) const {
    return if_eq_(_src, __Len);
  }
  bool operator==(const char* _src) const;
  bool operator==(const std::string& _src) const;
  bool operator==(const FortranStringView& _src) const;

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
   * @param _begin1  One-based index for the beginning index.
   * @param _back1   One-based index for the inclusive ending index.
   * @return         A new FortranStringView object.
   */
  FortranStringView operator()(int _begin1, int _back1) const;
};

typedef FortranStringView fstr_view;
TINKER_NAMESPACE_END

#endif
