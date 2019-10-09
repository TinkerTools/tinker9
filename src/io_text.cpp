#include "io_text.h"
#include <algorithm>
#include <bitset>
#include <cassert>
#include <regex>
#include <sstream>

TINKER_NAMESPACE_BEGIN
bool Text::is_ws(char ch)
{
   static std::string ws = whitespaces;
   return ws.find(ch) != std::string::npos;
}

void Text::replace(std::string& src, std::string old, char r)
{
   std::bitset<128> alphabets;
   for (char c : old) {
      alphabets[c] = 1;
   }
   for (char& rc : src) {
      if (alphabets[rc]) {
         rc = r;
      }
   }
}

void Text::replace_by_kv(std::string& src, std::string key, std::string value)
{
   std::regex rpl(key);
   src = std::regex_replace(src, rpl, value);
}

std::vector<std::string> Text::split(std::string str, std::string delimiters)
{
   assert(delimiters.size());

   char dlm = delimiters[0];
   replace(str, delimiters, dlm);
   std::istringstream ssin(str);
   std::vector<std::string> vs;
   std::string line;

   while (std::getline(ssin, line, dlm)) {
      if (line != "") {
         vs.push_back(line);
      }
   }

   return vs;
}

void Text::upcase(std::string& s)
{
   std::transform(s.begin(), s.end(), s.begin(), ::toupper);
}

void Text::lowcase(std::string& s)
{
   std::transform(s.begin(), s.end(), s.begin(), ::tolower);
}
TINKER_NAMESPACE_END
