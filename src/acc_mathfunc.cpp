#include "acc_mathfunc.h"
#include <cassert>

TINKER_NAMESPACE_BEGIN
namespace parallel {
template <class T>
T reduce_sum_tmpl(const T* gpu_a, size_t cpu_n) {
  T val = 0;
  #pragma acc parallel loop deviceptr(gpu_a) reduction(+:val)
  for (size_t i = 0; i < cpu_n; ++i)
    val += gpu_a[i];
  return val;
}

int reduce_sum(const int* a, size_t n) { return reduce_sum_tmpl(a, n); }

float reduce_sum(const float* a, size_t n) { return reduce_sum_tmpl(a, n); }

double reduce_sum(const double* a, size_t n) { return reduce_sum_tmpl(a, n); }

unsigned long long reduce_sum(const unsigned long long* a, size_t n) {
  return reduce_sum_tmpl(a, n);
}

//====================================================================//

template <class T>
void reduce_sum2_tmpl(T* __restrict__ h_ans, size_t hn, const T* __restrict__ v,
                      size_t nelem, size_t neach) {
  assert(hn <= neach);
  #pragma acc parallel loop gang deviceptr(v)
  for (size_t iv = 0; iv < hn; ++iv) {
    T ans = 0;
    #pragma acc loop vector reduction(+:ans)
    for (size_t ig = 0; ig < nelem; ++ig)
      ans += v[iv + ig * neach];
    h_ans[iv] = ans;
  }
}

void reduce_sum2(int* h_ans, size_t hn, const int* v, size_t nelem,
                 size_t neach) {
  reduce_sum2_tmpl(h_ans, hn, v, nelem, neach);
}

void reduce_sum2(float* h_ans, size_t hn, const float* v, size_t nelem,
                 size_t neach) {
  reduce_sum2_tmpl(h_ans, hn, v, nelem, neach);
}

void reduce_sum2(double* h_ans, size_t hn, const double* v, size_t nelem,
                 size_t neach) {
  reduce_sum2_tmpl(h_ans, hn, v, nelem, neach);
}

void reduce_sum2(unsigned long long* h_ans, size_t hn,
                 const unsigned long long* v, size_t nelem, size_t neach) {
  reduce_sum2_tmpl(h_ans, hn, v, nelem, neach);
}

//====================================================================//

template <class T>
T dotprod_acc_impl(const T* __restrict__ gpu_a, const T* __restrict__ gpu_b,
                   size_t cpu_n) {
  T val = 0;
  #pragma acc parallel loop deviceptr(gpu_a,gpu_b) reduction(+:val)
  for (size_t i = 0; i < cpu_n; ++i)
    val += gpu_a[i] * gpu_b[i];
  return val;
}

float dotprod(const float* a, const float* b, size_t n) {
  return dotprod_acc_impl(a, b, n);
}

double dotprod(const double* a, const double* b, size_t n) {
  return dotprod_acc_impl(a, b, n);
}

//====================================================================//

template <class T>
void scale_array_acc_impl(T* gpu_dst, T scal, size_t nelem) {
  #pragma acc parallel loop deviceptr(gpu_dst)
  for (size_t i = 0; i < nelem; ++i) {
    gpu_dst[i] *= scal;
  }
}

void scale_array(float* gpu_dst, float scal, size_t nelem) {
  scale_array_acc_impl(gpu_dst, scal, nelem);
}

void scale_array(double* gpu_dst, double scal, size_t nelem) {
  scale_array_acc_impl(gpu_dst, scal, nelem);
}
}
TINKER_NAMESPACE_END
