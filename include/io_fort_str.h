#ifndef TINKER_IO_FORT_STR_H_
#define TINKER_IO_FORT_STR_H_

#include "macro.h"
#include <string>

TINKER_NAMESPACE_BEGIN
class FortranStringView;
typedef FortranStringView fstr_view;

/**
 * @brief
 * references to the non-NULL-terminated fortran strings that provides a few
 * @c std::string -like methods to handle the fortran strings in C++ programs
 */
class FortranStringView {
private:
  char* const b_; ///< begin in [begin, end)
  char* const e_; ///< end in [begin, end)
  FortranStringView() = delete;

  /**
   * @brief
   * if @c dst != @c src, copy the @c first_n characters from @c src to @c dst;
   * fill @c dst with blanks if @c first_n is less than @c dstlen;
   * @c dst is NOT NULL-terminated NULL.
   */
  static void copy_with_blank_(char* dst, size_t dstlen, const char* src,
                               size_t first_n);

  /**
   * @brief
   * compare to string @c src of @c len;
   * the shorter string is filled by blanks prior to comparison
   */
  bool if_eq_(const char* src, size_t len) const;

public:
  template <size_t Len>
  FortranStringView(const char (&src)[Len])
      : b_(const_cast<char*>(src))
      , e_(b_ + Len) {}
  FortranStringView(const char* src, size_t len);
  FortranStringView(const char* src);
  FortranStringView(const std::string& src);

  // assignment
  template <size_t Len>
  FortranStringView& operator=(const char (&src)[Len]) {
    copy_with_blank_(b_, size(), src, Len);
    return *this;
  }
  FortranStringView& operator=(const char* src);
  FortranStringView& operator=(const std::string& src);
  FortranStringView& operator=(const FortranStringView& src);

  // comparison
  template <size_t Len>
  bool operator==(const char (&src)[Len]) const {
    return if_eq_(src, Len);
  }
  bool operator==(const char* src) const;
  bool operator==(const std::string& src) const;
  bool operator==(const FortranStringView& src) const;

  /// @return
  /// max number of characters in the string
  size_t size() const;

  /// @return
  /// length of string, ignoring any trailing blanks
  size_t len_trim() const;

  /// @return
  /// trimmed result in std::string.
  std::string trim() const;

  /**
   * @brief
   * analogous to fortran str(x:y) syntax
   *
   * @param begin1
   * one-based index for the beginning index
   *
   * @param back1
   * one-based index for the inclusive ending index
   *
   * @return
   * a new FortranStringView object
   */
  FortranStringView operator()(int begin1, int back1) const;

  /**
   * @brief
   * analogous to fortran str(x:) syntax
   *
   * @param begin1
   * one-based index for the beginning index
   *
   * @return
   * a new FortranStringView object
   */
  FortranStringView operator()(int begin1) const;
};
TINKER_NAMESPACE_END

#endif
