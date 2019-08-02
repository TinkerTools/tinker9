#include "acc_seq.h"

TINKER_NAMESPACE_BEGIN
// TODO: update dotprod implementation
template <class T>
T dotprod_acc_impl(const T* gpu_a, const T* gpu_b, int cpu_n) {
  T val = 0;
  #pragma acc parallel loop independent copy(val) reduction(+:val)\
              deviceptr(gpu_a,gpu_b)
  for (int i = 0; i < cpu_n; ++i) {
    val += gpu_a[i] * gpu_b[i];
  }
  return val;
}

float dotprod(const float* a, const float* b, int n) {
  return dotprod_acc_impl<float>(a, b, n);
}

double dotprod(const double* a, const double* b, int n) {
  return dotprod_acc_impl<double>(a, b, n);
}

template <class T>
void scale_array_acc_impl(T* gpu_dst, T scal, int nelem) {
  #pragma acc parallel loop independent deviceptr(gpu_dst)
  for (int i = 0; i < nelem; ++i) {
    gpu_dst[i] *= scal;
  }
}

void scale_array(float* gpu_dst, float scal, int nelem) {
  scale_array_acc_impl<float>(gpu_dst, scal, nelem);
}

void scale_array(double* gpu_dst, double scal, int nelem) {
  scale_array_acc_impl<double>(gpu_dst, scal, nelem);
}
TINKER_NAMESPACE_END
