#pragma once
#include <type_traits>

namespace tinker {
/// \ingroup cuda_syntax
/// \brief Used as `eq<T1,T2>()` for two type identifiers.
template <class T, class U>
__device__
__host__
constexpr bool eq()
{
   return std::is_same<T, U>::value;
}
}
