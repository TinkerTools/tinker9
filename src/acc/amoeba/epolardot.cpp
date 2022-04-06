#include "ff/amoebamod.h"
#include "ff/atom.h"
#include "ff/elec.h"
#include "seq/add.h"

namespace tinker {
void epolar0DotProd_acc(const real (*gpu_uind)[3], const real (*gpu_udirp)[3])
{
   const real f = -0.5 * electric / dielec;

   auto bufsize = bufferSize();

   #pragma acc parallel loop independent async\
               deviceptr(ep,gpu_uind,gpu_udirp,polarity_inv)
   for (int i = 0; i < n; ++i) {
      int offset = i & (bufsize - 1);
      real e = polarity_inv[i] *
         (gpu_uind[i][0] * gpu_udirp[i][0] + gpu_uind[i][1] * gpu_udirp[i][1] +
            gpu_uind[i][2] * gpu_udirp[i][2]);
      atomic_add(f * e, ep, offset);
   }
}
}
