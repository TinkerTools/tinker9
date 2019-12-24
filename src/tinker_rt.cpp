#include "tinker_rt.h"
#include <algorithm>
#include <tinker/detail/argue.hh>

TINKER_NAMESPACE_BEGIN
void nextarg(size_t len, char* str, int& exist)
{
   const char blank = ' ';
   std::memset(str, blank, len);
   exist = false;

   if (argue::narg != 0) {
      size_t length = std::min(len, sizeof(argue::arg[argue::maxarg]));
      for (size_t i = 1; i <= argue::narg; ++i) {
         if (argue::listarg[i]) {
            argue::listarg[i] = false;
            std::strncpy(str, argue::arg[i], length);
            exist = true;
            break;
         }
      }
   }
}


namespace detail {
template <class T>
void vstr_to_val(T&, std::string);


template <>
void vstr_to_val<int>(int& v, std::string vstr)
{
   v = std::stoi(vstr);
}
}


template <class T>
void get_kv_pair(std::string k, T& v, T vdefault)
{
   std::string value_str = "";
   for (int i = 0; i < keys::nkey; ++i) {
      fstr_view record = keys::keyline[i];
      auto vs = Text::split(record.trim());
      if (vs.size()) {
         std::string keyword = vs.at(0);
         Text::upcase(keyword);
         if (keyword == k && vs.size() > 1) {
            value_str = vs[1];
         }
      }
   }
   if (value_str != "") {
      detail::vstr_to_val(v, value_str);
   } else {
      v = vdefault;
   }
}
template void get_kv_pair<int>(std::string, int&, int);
TINKER_NAMESPACE_END
