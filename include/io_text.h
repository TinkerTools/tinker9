#pragma once
#include "macro.h"
#include <string>
#include <vector>


namespace tinker {
/**
 * \ingroup io
 * \brief Plain ascii text stored by lines.
 */
class Text : public std::vector<std::string>
{
public:
   /// \brief White space characters.
   static constexpr const char* whitespaces = " \t\n\v\f\r";


   /// \return True if `ch` is one of the white space characters.
   static bool is_ws(char ch);


   /// \return An `std::string` object from a c-style `char` array.
   template <size_t Len>
   static std::string string(const char (&src)[Len])
   {
      return std::string(&src[0], &src[0] + Len);
   }


   /// \brief In string `src`, replace any characters appeared in `old` to `r`.
   static void replace(std::string& src, std::string old, char r);


   /// \brief In string `src`, replace any string `key` to string `value`.
   static void replace_by_kv(std::string& src, std::string key,
                             std::string value);


   /// \brief
   /// Split a string `str` by `delimiters` and return a vector of strings.
   static std::vector<std::string> split(std::string str,
                                         std::string delimiters = whitespaces);


   /// \brief Transform a string to upper case.
   static void upcase(std::string&);


   /// \brief Transform a string to lower case.
   static void lowcase(std::string&);
};
}
