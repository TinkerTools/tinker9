#pragma once
#include "macro.h"
#include <type_traits>


namespace tinker {
template <class T, class U>
__device__
__host__
constexpr bool eq()
{
   return std::is_same<T, U>::value;
}
}
