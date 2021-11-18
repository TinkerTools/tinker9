#pragma once
#include "macro.h"
#include <string>


namespace tinker {
/**
 * \ingroup io
 * Reference to the non-NULL-terminated fortran strings,
 * which provides a few `std::string`-style methods to handle
 * the fortran strings in C++ programs.
 */
class FortranStringView
{
private:
   char* const m_b; // begin in [begin, end)
   char* const m_e; // end in [begin, end)


   FortranStringView() = delete; // no default constructor


   /**
    * If `dst != src`, copy the `first_n` characters from `src` to `dst`;
    * fill `dst` with blanks if `first_n` is less than `dstlen`;
    * `dst` is NOT NULL-terminated.
    */
   static void copy_with_blank(char* dst, size_t dstlen, const char* src,
                               size_t first_n);


   /**
    * Compare to a string `src` of `len`.
    * The shorter string is filled by blanks prior to comparison.
    */
   bool if_eq(const char* src, size_t len) const;


   /**
    * Returns the max number of characters in the string including the trailing
    * whitespaces.
    */
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


   /**
    * Analogous to Fortran str(x:y) syntax.
    * \param begin1  One-based index for the beginning index.
    * \param back1   One-based index for the inclusive ending index.
    * \return        A new FortranStringView object.
    */
   FortranStringView operator()(int begin1, int back1) const;


   /**
    * Analogous to Fortran str(x:) syntax.
    * \param begin1  One-based index for the beginning index.
    * \return        A new FortranStringView object
    */
   FortranStringView operator()(int begin1) const;
};


/**
 * \typedef fstr_view
 * \ingroup io
 * An type alias of FortranStringView.
 */
using fstr_view = FortranStringView;
}
