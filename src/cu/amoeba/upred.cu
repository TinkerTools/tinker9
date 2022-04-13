#include "seq/launch.h"

namespace tinker {
__global__
void ulspredSaveP1_cu1(int n, real (*restrict ud)[3], real (*restrict up)[3],
   const real (*restrict uind)[3], const real (*restrict uinp)[3])
{
   if (uinp) {
      for (int i = ITHREAD; i < n; i += STRIDE) {
         ud[i][0] = uind[i][0];
         ud[i][1] = uind[i][1];
         ud[i][2] = uind[i][2];
         up[i][0] = uinp[i][0];
         up[i][1] = uinp[i][1];
         up[i][2] = uinp[i][2];
      }
   } else {
      for (int i = ITHREAD; i < n; i += STRIDE) {
         ud[i][0] = uind[i][0];
         ud[i][1] = uind[i][1];
         ud[i][2] = uind[i][2];
      }
   }
}

void ulspredSaveP1_cu(real (*ud)[3], real (*up)[3], const real (*uind)[3], const real (*uinp)[3])
{
   launch_k1s(g::s0, n, ulspredSaveP1_cu1, n, ud, up, uind, uinp);
}
}
