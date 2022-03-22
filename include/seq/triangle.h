#pragma once
#include "math/inc.h"

namespace tinker {
__device__
inline int xy_to_tri(int x, int y)
{
   long long lx = x;
   int base = (lx + 1) * lx / 2;
   return base + y;
}

__device__
inline void tri_to_xy(int f, int& x, int& y)
{
   long long lf = 8ll * f + 1;
   double ff = lf;
   double fa = (sqrt(ff) - 1) / 2;
   x = fa;
   y = f - xy_to_tri(x, 0);
}
}
