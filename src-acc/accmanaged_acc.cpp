#include "accmanaged.h"


namespace tinker {
void accmanaged_data(rc_op op)
{
   if (op & rc_dealloc) {
      // #pragma acc exit data async delete()
   }


   if (op & rc_alloc) {
      // #pragma acc enter data async create()
   }
}
}
