#include "tinker_rt.h"
#include <algorithm>
#include <tinker/detail/argue.hh>
#include <tinker/detail/keys.hh>

namespace tinker {
void nextarg(size_t len, char* str, int& exist)
{
   const char blank = ' ';
   std::memset(str, blank, len);
   exist = false;

   if (argue::narg != 0) {
      size_t length = std::min(len, sizeof(argue::arg[argue::maxarg]));
      for (int i = 1; i <= argue::narg; ++i) {
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

template <>
void vstr_to_val<double>(double& v, std::string vstr)
{
   v = std::stod(vstr);
}

template <>
void vstr_to_val<std::string>(std::string& v, std::string vstr)
{
   v = vstr;
   Text::upcase(v);
}
}

template <class T1, class T2>
void get_kv(std::string k, T1& v, T2 vdefault)
{
   std::string value_str = "";
   for (int i = 0; i < keys::nkey; ++i) {
      FstrView record = keys::keyline[i];
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
template void get_kv(std::string, int&, int);
template void get_kv(std::string, double&, double);
template void get_kv(std::string, std::string&, std::string);
template void get_kv(std::string, std::string&, const char*);
template <>
void get_kv<std::vector<std::string>, const char*>(
   std::string k, std::vector<std::string>& v, const char*)
{
   for (int i = 0; i < keys::nkey; ++i) {
      FstrView record = keys::keyline[i];
      auto vs = Text::split(record.trim());
      if (vs.size()) {
         std::string keyword = vs.at(0);
         Text::upcase(keyword);
         if (keyword == k && vs.size() > 1) {
            v.clear();
            for (size_t j = 1; j < vs.size(); ++j) {
               v.push_back(vs[j]);
               Text::upcase(v[j - 1]);
            }
         }
      }
   }
}

template <class T>
void get_kbool(std::string k, T& v, bool v_if_k_not_found)
{
   bool found = false;
   for (int i = 0; i < keys::nkey; ++i) {
      FstrView record = keys::keyline[i];
      auto vs = Text::split(record.trim());
      if (vs.size()) {
         std::string keyword = vs.at(0);
         Text::upcase(keyword);
         if (keyword == k) {
            found = true;
         }
      }
   }
   if (found)
      v = !v_if_k_not_found;
   else
      v = v_if_k_not_found;
}
template void get_kbool(std::string k, bool& v, bool vdefault);
template void get_kbool(std::string k, int& v, bool vdefault);
}
