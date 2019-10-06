#ifndef TINKER_IO_TEXT_H_
#define TINKER_IO_TEXT_H_

#include "macro.h"
#include <string>
#include <vector>

TINKER_NAMESPACE_BEGIN
/// @brief
/// plain ascii text stored by lines
class Text : public std::vector<std::string>
{
public:
   /// @brief
   /// white space characters
   static constexpr const char* whitespaces = " \t\n\v\f\r";

   /// @return
   /// @c true if @c ch is one of the white space characters
   static bool is_ws (char ch);

   /// @return
   /// an @c std::string object from a c-style @c char array
   template <size_t Len>
   static std::string string (const char (&src)[Len])
   {
      return std::string (&src[0], &src[0] + Len);
   }

   /// @brief
   /// in string @c src, replace any characters appeared in @c old to char @c r
   static void replace (std::string& src, std::string old, char r);

   /// @brief
   /// in string @c src, replace any string @c key to string @c value
   static void replace_by_kv (std::string& src, std::string key,
                              std::string value);

   /// @brief
   /// split a string @c str to a vector of strings by @c delimiters
   static std::vector<std::string> split (std::string str,
                                          std::string delimiters = whitespaces);

   /// @brief
   /// transform a string to upper or lower case
   /// @{
   static void upcase (std::string&);
   static void lowcase (std::string&);
   /// @}
};
TINKER_NAMESPACE_END

#endif
