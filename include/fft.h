#pragma once
#include "macro.h"


TINKER_NAMESPACE_BEGIN
struct FFTPlan
{
   template <class T>
   T& self()
   {
      return *reinterpret_cast<T*>(this);
   }
};
TINKER_NAMESPACE_END
