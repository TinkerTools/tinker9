#pragma once
#include "macro.h"
#include <type_traits>


namespace tinker {
/// \ingroup acc_syntax
template <class T, class U>
constexpr bool eq()
{
   return std::is_same<T, U>::value;
}
}
