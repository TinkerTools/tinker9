#pragma once
#include "macro.h"


TINKER_NAMESPACE_BEGIN
/// \ingroup math
/// \see max_of(T, T2, Ts...)
template <class T, class T2>
T max_of(T a, T2 b)
{
   return (a < b) ? b : a;
}


/// \ingroup math
/// \return
/// The maximum value from a variadic list of numbers.
template <class T, class T2, class... Ts>
T max_of(T a, T2 b, Ts... cs)
{
   return max_of(max_of(a, b), cs...);
}


/// \ingroup math
/// \see min_of(T, T2, Ts...)
template <class T, class T2>
T min_of(T a, T2 b)
{
   return (a < b) ? a : b;
}


/// \ingroup math
/// \return
/// The minimum value from a variadic list of numbers.
template <class T, class T2, class... Ts>
T min_of(T a, T2 b, Ts... cs)
{
   return min_of(min_of(a, b), cs...);
}
TINKER_NAMESPACE_END
