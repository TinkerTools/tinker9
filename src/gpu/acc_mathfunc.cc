#include "gpu/acc_mathfunc.h"

TINKER_NAMESPACE_BEGIN
template <class T>
void dotprod_acc_impl(T* cpu_ans, const T* gpu_a, const T* gpu_b, int cpu_n) {
  T val = 0;
  #pragma acc parallel loop independent copy(val) reduction(+:val)\
              deviceptr(gpu_a,gpu_b)
  for (int i = 0; i < cpu_n; ++i) {
    val += gpu_a[i] * gpu_b[i];
  }
  *cpu_ans = val;
}

void dotprod(float* ans, const float* a, const float* b, int n) {
  dotprod_acc_impl<float>(ans, a, b, n);
}

void dotprod(double* ans, const double* a, const double* b, int n) {
  dotprod_acc_impl<double>(ans, a, b, n);
}
TINKER_NAMESPACE_END
