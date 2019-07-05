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
   * If dst_ != src_, copy the first_n_ characters from src_ to dst_;
   * fill dst_ with blanks if first_n_ is less than dstlen_;
   * dst_ is NOT terminated by NULL.
   */
  static void copy_with_blank_(char* dst_, size_t dstlen_, const char* src_,
                               size_t first_n_);

  /**
   * Compare to string src_ of len_.
   * The shorter string is filled by blanks prior to comparison.
   */
  bool if_eq_(const char* src_, size_t len_) const;

public:
  //
  // Constructors.
  //
  template <size_t Len_>
  FortranStringView(const char (&src_)[Len_])
      : b_(const_cast<char*>(src_)), e_(b_ + Len_) {}
  FortranStringView(const char* src_, size_t len_);
  FortranStringView(const char* src_);
  FortranStringView(const std::string& src_);

  //
  // Assignment operators.
  //
  template <size_t Len_>
  FortranStringView& operator=(const char (&src_)[Len_]) {
    copy_with_blank_(b_, size(), src_, Len_);
    return *this;
  }
  FortranStringView& operator=(const char* src_);
  FortranStringView& operator=(const std::string& src_);
  FortranStringView& operator=(const FortranStringView& src_);

  //
  // Comparisons.
  //
  template <size_t Len_>
  bool operator==(const char (&src_)[Len_]) const {
    return if_eq_(src_, Len_);
  }
  bool operator==(const char* src_) const;
  bool operator==(const std::string& src_) const;
  bool operator==(const FortranStringView& src_) const;

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
   * @param begin1_  One-based index for the beginning index.
   * @param back1_   One-based index for the inclusive ending index.
   * @return         A new FortranStringView object.
   */
  FortranStringView operator()(int begin1_, int back1_) const;
};

typedef FortranStringView fstr_view;
TINKER_NAMESPACE_END

#endif
