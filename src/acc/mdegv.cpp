#include "mdegv.h"
#include "add.h"
#include "mdpq.h"

namespace tinker {
void scale_gradient_acc(double scale, grad_prec* g0x, grad_prec* g0y, grad_prec* g0z)
{
   real s = scale;
#if TINKER_DETERMINISTIC_FORCE
   #pragma acc parallel loop independent async\
           deviceptr(g0x,g0y,g0z)
   for (int i = 0; i < n; ++i) {
      real dx = s * to_flt_acc<real>(g0x[i]);
      real dy = s * to_flt_acc<real>(g0y[i]);
      real dz = s * to_flt_acc<real>(g0z[i]);
      g0x[i] = cvt_to<grad_prec>(dx);
      g0y[i] = cvt_to<grad_prec>(dy);
      g0z[i] = cvt_to<grad_prec>(dz);
   }
#else
   #pragma acc parallel loop independent async\
           deviceptr(g0x,g0y,g0z)
   for (int i = 0; i < n; ++i) {
      g0x[i] *= s;
      g0y[i] *= s;
      g0z[i] *= s;
   }
#endif
}

void sum_gradient_acc(grad_prec* g0x, grad_prec* g0y, grad_prec* g0z, const grad_prec* g1x,
   const grad_prec* g1y, const grad_prec* g1z)
{
   #pragma acc parallel loop independent async\
           deviceptr(g0x,g0y,g0z,g1x,g1y,g1z)
   for (int i = 0; i < n; ++i) {
      g0x[i] += g1x[i];
      g0y[i] += g1y[i];
      g0z[i] += g1z[i];
   }
}

void sum_gradient_acc(double ss, grad_prec* g0x, grad_prec* g0y, grad_prec* g0z,
   const grad_prec* g1x, const grad_prec* g1y, const grad_prec* g1z)
{
   real s = ss;
#if TINKER_DETERMINISTIC_FORCE
   #pragma acc parallel loop independent async\
           deviceptr(g0x,g0y,g0z,g1x,g1y,g1z)
   for (int i = 0; i < n; ++i) {
      real dx = s * to_flt_acc<real>(g1x[i]);
      real dy = s * to_flt_acc<real>(g1y[i]);
      real dz = s * to_flt_acc<real>(g1z[i]);
      g0x[i] += cvt_to<grad_prec>(dx);
      g0y[i] += cvt_to<grad_prec>(dy);
      g0z[i] += cvt_to<grad_prec>(dz);
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
}
