#include "tool/iotext.h"
#include <algorithm>
#include <bitset>
#include <cassert>
#include <sstream>

namespace tinker {
bool Text::isWhiteSpace(char ch)
{
   static std::string ws = whitespaces;
   return ws.find(ch) != std::string::npos;
}

void Text::replace(std::string& src, std::string old, char r)
{
   std::bitset<256> alphabets;
   for (char ic : old) {
      unsigned char c = ic;
      alphabets[c] = true;
   }
   for (char& rc : src) {
      unsigned char c = rc;
      if (alphabets[c]) {
         rc = r;
      }
   }
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
}
