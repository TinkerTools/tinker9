#pragma once
#include "macro.h"

namespace tinker {
/// \ingroup math
/// \see maxOf(T, T2, Ts...)
template <class T, class T2>
T maxOf(T a, T2 b)
{
   return (a < b) ? b : a;
}

/// \ingroup math
/// \return The maximum value from a variadic list of numbers.
template <class T, class T2, class... Ts>
T maxOf(T a, T2 b, Ts... cs)
{
   return maxOf(maxOf(a, b), cs...);
}

/// \ingroup math
/// \see minOf(T, T2, Ts...)
template <class T, class T2>
T minOf(T a, T2 b)
{
   return (a < b) ? a : b;
}

/// \ingroup math
/// \return The minimum value from a variadic list of numbers.
template <class T, class T2, class... Ts>
T minOf(T a, T2 b, Ts... cs)
{
   return minOf(minOf(a, b), cs...);
}
}
