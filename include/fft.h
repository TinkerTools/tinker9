#pragma once
#include "macro.h"
#include <type_traits>


namespace tinker {
struct FFTPlan
{
   template <class T>
   T& self()
   {
      static_assert(std::is_base_of<FFTPlan, T>::value, "");
      return *reinterpret_cast<T*>(this);
   }

   virtual ~FFTPlan() {}
};
}
