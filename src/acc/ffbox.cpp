#include "ff/box.h"

namespace tinker {
void boxData_acc(RcOp op)
{
   if (op & rc_dealloc) {
      #pragma acc exit data async delete(lvec1,lvec2,lvec3,recipa,recipb,recipc)
   }

   if (op & rc_alloc) {
      #pragma acc enter data async create(lvec1,lvec2,lvec3,recipa,recipb,recipc)
   }
}

void boxCopyin_acc()
{
   #pragma acc update async device(lvec1,lvec2,lvec3,recipa,recipb,recipc)
}
}
