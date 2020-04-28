#include "add.h"
#include "mdegv.h"
#include "mdpq.h"
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


void sum_gradient_acc(grad_prec* g0x, grad_prec* g0y, grad_prec* g0z,
                      double scale, const grad_prec* g1x, const grad_prec* g1y,
                      const grad_prec* g1z)
{
   real s = scale;
#if TINKER_DETERMINISTIC_FORCE
   #pragma acc parallel loop independent async\
           deviceptr(g0x,g0y,g0z,g1x,g1y,g1z)
   for (int i = 0; i < n; ++i) {
      real dx = s * to_flt_acc<real>(g1x[i]);
      real dy = s * to_flt_acc<real>(g1y[i]);
      real dz = s * to_flt_acc<real>(g1z[i]);
      g0x[i] += acc_to<grad_prec>(dx);
      g0y[i] += acc_to<grad_prec>(dy);
      g0z[i] += acc_to<grad_prec>(dz);
   }
#else
   #pragma acc parallel loop independent async\
           deviceptr(g0x,g0y,g0z,g1x,g1y,g1z)
   for (int i = 0; i < n; ++i) {
      g0x[i] += s * g1x[i];
      g0y[i] += s * g1y[i];
      g0z[i] += s * g1z[i];
   }
#endif
}
TINKER_NAMESPACE_END
