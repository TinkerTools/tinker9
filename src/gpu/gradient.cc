#include "gpu/data.h"
#include <cstdio>

TINKER_NAMESPACE_BEGIN
namespace gpu {
static const int N = 16000000;
static real x[N];

extern "C" {
void tinker_gpu_gradient() {
  for (int i = N - 100; i < N; ++i) {
    x[i] = -i;
  }

#pragma acc kernels
  {
    for (int i = 0; i < N; ++i) {
      x[i] += (i + 100);
    }
  }

  for (int i = N - 100; i < N; ++i) {
    printf("index = %d, val = %lf\n", i, x[i]);
  }
}
}
}
TINKER_NAMESPACE_END
