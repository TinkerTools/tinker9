#include "tool/gpucard.h"

namespace tinker {
void gpuData(RcOp op)
{
   if (op & rc_dealloc) {
      ndevice = 0;
      idevice = -1;
   }

   if (op & rc_init) {
      ndevice = 1;
      idevice = 0;
   }
}

int gpuGridSize(int)
{
   return 1;
}

int gpuMaxNParallel(int)
{
   return 1;
}
}
