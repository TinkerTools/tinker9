#pragma once
#include "macro.h"
#include <string>
#include <vector>


namespace tinker {
/**
 * \ingroup io
 * Plain ascii text stored by lines.
 */
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
   static void replace_by_kv(std::string& src, std::string key,
                             std::string value);


   /// Splits a string `str` by `delimiters` and returns a vector of strings.
   static std::vector<std::string> split(std::string str,
                                         std::string delimiters = whitespaces);


   /// Transforms a string to upper case.
   static void upcase(std::string&);


   /// Transforms a string to lower case.
   static void lowcase(std::string&);
};
}
