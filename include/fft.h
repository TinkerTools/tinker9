#pragma once
#include "macro.h"
#include <type_traits>


TINKER_NAMESPACE_BEGIN
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
TINKER_NAMESPACE_END
