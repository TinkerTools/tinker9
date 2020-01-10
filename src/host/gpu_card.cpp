#include "gpu_card.h"


TINKER_NAMESPACE_BEGIN
void gpu_card_data(rc_op op)
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


int get_grid_size(int)
{
   return 1;
}
TINKER_NAMESPACE_END
