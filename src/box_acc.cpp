#include "box.h"


namespace tinker {
void box_data_acc(rc_op op)
{
   if (op & rc_dealloc) {
      #pragma acc exit data delete(lvec1,lvec2,lvec3,recipa,recipb,recipc)
   }


   if (op & rc_alloc) {
      #pragma acc enter data create(lvec1,lvec2,lvec3,recipa,recipb,recipc)
   }
}


void box_copyin_acc()
{
   #pragma acc update device(lvec1,lvec2,lvec3,recipa,recipb,recipc)
}
}
