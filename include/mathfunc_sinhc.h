#pragma once
#include <cmath>

namespace tinker {
/**
 * \ingroup math
 */
template <class T>
T sinhc(T x)
{
   T y = std::fabs(x);
   if (y <= (T)1.0e-8)
      return (T)1;
   else
      return std::sinh(y) / y;
}
}
