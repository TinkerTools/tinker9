#include "tool/gpu_card.h"


namespace tinker {
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


int gpu_max_nparallel(int)
{
   return 1;
}
}
