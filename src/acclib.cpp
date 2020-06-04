#include "tool/acclib.h"


namespace tinker {
void wait_queue(LPFlag flag)
{
   if ((flag & LPFlag::WAIT) && !(flag & LPFlag::DEFAULT_Q)) {
      #pragma acc wait(async_queue)
   }
}
}
