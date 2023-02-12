#include "seq/rotpole.h"
#include "ff/modamoeba.h"
#include "ff/atom.h"
#include "math/libfunc.h"

namespace tinker {
void chkpole_acc()
{
   #pragma acc parallel loop independent async deviceptr(x,y,z,zaxis,pole)
   for (int i = 0; i < n; ++i)
      chkpoleAtomI(i, pole, zaxis, x, y, z);
}

void rotpole_acc()
{
   #pragma acc parallel loop independent async deviceptr(x,y,z,zaxis,rpole,pole)
   for (int i = 0; i < n; ++i)
      rotpoleAtomI(i, rpole, pole, zaxis, x, y, z);
}

void rotrepole_acc()
{
   #pragma acc parallel loop independent async deviceptr(x,y,z,zaxis,rrepole,repole)
   for (int i = 0; i < n; ++i)
      rotpoleAtomI(i, rrepole, repole, zaxis, x, y, z);
}
}
