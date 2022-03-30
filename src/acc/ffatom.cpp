#include "ff/atom.h"

namespace tinker {
void copyPosToXyz_acc()
{
   if CONSTEXPR (sizeof(pos_prec) == sizeof(real))
      return;

   #pragma acc parallel loop independent async deviceptr(x,y,z,xpos,ypos,zpos)
   for (int i = 0; i < n; ++i) {
      x[i] = xpos[i];
      y[i] = ypos[i];
      z[i] = zpos[i];
   }
}
}
