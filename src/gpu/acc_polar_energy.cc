#include "gpu/e_polar.h"
#include "mod_md.h"

TINKER_NAMESPACE_BEGIN
void epolar0_dotprod(const real (*gpu_uind)[3], const real (*gpu_udirp)[3]) {
  real e = 0;
  #pragma acc parallel loop independent copy(e) reduction(+:e)\
              deviceptr(gpu_uind,gpu_udirp,polarity_inv)
  for (int i = 0; i < n; ++i) {
    e += polarity_inv[i] *
        (gpu_uind[i][0] * gpu_udirp[i][0] + gpu_uind[i][1] * gpu_udirp[i][1] +
         gpu_uind[i][2] * gpu_udirp[i][2]);
  }

  const real f = -0.5 * electric / dielec;
  #pragma acc serial deviceptr(ep)
  { *ep = f * e; }
}
TINKER_NAMESPACE_END
