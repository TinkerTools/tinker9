#include "tool/gpucard.h"

namespace tinker {
void gpuData(RcOp op)
{
   if (op & RcOp::DEALLOC) {
      ndevice = 0;
      idevice = -1;
   }

   if (op & RcOp::INIT) {
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
