#include "accmanaged.h"


namespace tinker {
void accmanaged_data(rc_op op)
{
   using namespace detail;


   if (op & rc_dealloc) {
      // #pragma acc exit data async delete()
   }


   if (op & rc_alloc) {
      // #pragma acc enter data async create()
   }
}


namespace detail {}
}
