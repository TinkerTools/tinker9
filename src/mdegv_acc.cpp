#include "mdegv.h"
#include "wait_queue.h"


TINKER_NAMESPACE_BEGIN
void zero_gradient_acc(DMFlag flag, size_t nelem, real* gx, real* gy, real* gz)
{
   bool sync = flag & DMFlag::DEFAULT_Q;
   if (sync) {
      #pragma acc parallel loop independent deviceptr(gx,gy,gz)
      for (size_t i = 0; i < nelem; ++i) {
         gx[i] = 0;
         gy[i] = 0;
         gz[i] = 0;
      }
   } else {
      #pragma acc parallel loop independent async deviceptr(gx,gy,gz)
      for (size_t i = 0; i < nelem; ++i) {
         gx[i] = 0;
         gy[i] = 0;
         gz[i] = 0;
      }
   }
   // if (flag & DMFlag::WAIT) {
   wait_queue(flag);
   // }
}


void zero_gradient_acc(DMFlag flag, size_t nelem, fixed* gx, fixed* gy,
                       fixed* gz)
{
   bool sync = flag & DMFlag::DEFAULT_Q;
   if (sync) {
      #pragma acc parallel loop independent deviceptr(gx,gy,gz)
      for (size_t i = 0; i < nelem; ++i) {
         gx[i] = 0;
         gy[i] = 0;
         gz[i] = 0;
      }
   } else {
      #pragma acc parallel loop independent async deviceptr(gx,gy,gz)
      for (size_t i = 0; i < nelem; ++i) {
         gx[i] = 0;
         gy[i] = 0;
         gz[i] = 0;
      }
   }
   // if (flag & DMFlag::WAIT) {
   wait_queue(flag);
   // }
}
TINKER_NAMESPACE_END
