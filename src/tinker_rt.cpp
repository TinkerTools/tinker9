#include "tinker_rt.h"
#include <algorithm>
#include <tinker/detail/argue.hh>


// For the hidden string length argument:
// GNU Fortran
// https://gcc.gnu.org/onlinedocs/gfortran/Argument-passing-conventions.html
// Intel Fortran
// https://software.intel.com/en-us/forums/intel-fortran-compiler/topic/637855
#if defined(TINKER_GFORTRAN) && (__GNUC__ <= 7)
extern "C" void TINKER_RT(evcorr1)(char*, double*, double*, int);
#else
extern "C" void TINKER_RT(evcorr1)(char*, double*, double*, size_t);
#endif


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


// evcorr.f
void evcorr1(char* mode, double* elrc, double* vlrc, int modelen)
{
   TINKER_RT(evcorr1)(mode, elrc, vlrc, modelen);
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
void get_kv_pair(std::string k, T1& v, T2 vdefault)
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
template void get_kv_pair(std::string, int&, int);
template void get_kv_pair(std::string, double&, double);
template void get_kv_pair(std::string, std::string&, std::string);
template void get_kv_pair(std::string, std::string&, const char*);
TINKER_NAMESPACE_END
