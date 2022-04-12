#include "seq/launch.h"

namespace tinker {
__global__
void copyPosToXyz_cu1(int n, const pos_prec* restrict xpoz, const pos_prec* restrict ypoz,
   const pos_prec* restrict zpoz, real* restrict x1, real* restrict y1, real* restrict z1)
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      x1[i] = xpoz[i];
      y1[i] = ypoz[i];
      z1[i] = zpoz[i];
   }
}

void copyPosToXyz_cu()
{
   launch_k1s(g::s0, n, copyPosToXyz_cu1, n, xpos, ypos, zpos, x, y, z);
}
}
