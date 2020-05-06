#pragma once
#include "macro.h"
#include <string>


namespace tinker {
class FortranStringView;
using fstr_view = FortranStringView;


/**
 * \ingroup io
 * \brief Reference to the non-NULL-terminated fortran strings,
 * which provides a few `std::string`-style methods to handle
 * the fortran strings in C++ programs.
 */
class FortranStringView
{
private:
   char* const b_; ///< begin in [begin, end)
   char* const e_; ///< end in [begin, end)
   FortranStringView() = delete;


   /**
    * \brief
    * If `dst != src`, copy the `first_n` characters from `src` to `dst`;
    * fill `dst` with blanks if `first_n` is less than `dstlen`;
    * `dst` is NOT NULL-terminated.
    */
   static void copy_with_blank_(char* dst, size_t dstlen, const char* src,
                                size_t first_n);


   /**
    * \brief Compare to a string `src` of `len`.
    * The shorter string is filled by blanks prior to comparison
    */
   bool if_eq_(const char* src, size_t len) const;


public:
   template <size_t Len>
   FortranStringView(const char (&src)[Len])
      : b_(const_cast<char*>(src))
      , e_(b_ + Len)
   {}
   FortranStringView(const char* src, size_t len);
   FortranStringView(const char* src);
   FortranStringView(const std::string& src);


   // assignment
   template <size_t Len>
   FortranStringView& operator=(const char (&src)[Len])
   {
      copy_with_blank_(b_, size(), src, Len);
      return *this;
   }
   FortranStringView& operator=(const char* src);
   FortranStringView& operator=(const std::string& src);
   FortranStringView& operator=(const FortranStringView& src);


   // comparison
   template <size_t Len>
   bool operator==(const char (&src)[Len]) const
   {
      return if_eq_(src, Len);
   }
   bool operator==(const char* src) const;
   bool operator==(const std::string& src) const;
   bool operator==(const FortranStringView& src) const;


   /// \return Max number of characters in the string.
   size_t size() const;


   /// \return Length of string, ignoring the trailing blanks.
   size_t len_trim() const;


   /// \return A string ignoring the trailing blanks.
   std::string trim() const;


   /**
    * \brief Analogous to Fortran str(x:y) syntax.
    * \param begin1 One-based index for the beginning index.
    * \param back1  One-based index for the inclusive ending index.
    * \return       A new FortranStringView object.
    */
   FortranStringView operator()(int begin1, int back1) const;


   /**
    * \brief Analogous to Fortran str(x:) syntax.
    * \param begin1 One-based index for the beginning index.
    * \return       A new FortranStringView object
    */
   FortranStringView operator()(int begin1) const;
};
}
