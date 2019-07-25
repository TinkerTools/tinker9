#include "util_text.h"
#include <algorithm>
#include <bitset>
#include <regex>
#include <sstream>

TINKER_NAMESPACE_BEGIN
void Text::replace(std::string& s, std::string old, char r) {
  std::bitset<128> alphabets;
  for (char c : old) {
    alphabets[c] = 1;
  }
  for (char& rc : s) {
    if (alphabets[rc]) {
      rc = r;
    }
  }
}

void Text::replace_by_kv(std::string& src, std::string key, std::string value) {
  std::regex rpl(key);
  src = std::regex_replace(src, rpl, value);
}

std::vector<std::string> Text::split(std::string str, std::string delimiters) {
  assert(delimiters.size());

  char delimiter = delimiters[0];
  replace(str, delimiters, delimiter);
  std::istringstream ssin(str);
  std::vector<std::string> vs;
  std::string line;

  while (std::getline(ssin, line, delimiter)) {
    if (line != "") {
      vs.push_back(line);
    }
  }

  return vs;
}

void Text::upcase(std::string& s) {
  std::transform(s.begin(), s.end(), s.begin(), ::toupper);
}

void Text::lowcase(std::string& s) {
  std::transform(s.begin(), s.end(), s.begin(), ::tolower);
}
TINKER_NAMESPACE_END
