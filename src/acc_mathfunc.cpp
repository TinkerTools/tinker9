#include "acc_mathfunc.h"
#include <cassert>

TINKER_NAMESPACE_BEGIN
template <class T>
T reduce_sum_tmpl(const T* gpu_a, int cpu_n) {
  T val = 0;
  #pragma acc parallel loop deviceptr(gpu_a) reduction(+:val)
  for (int i = 0; i < cpu_n; ++i)
    val += gpu_a[i];
  return val;
}

int reduce_sum(const int* a, int n) { return reduce_sum_tmpl(a, n); }

float reduce_sum(const float* a, int n) { return reduce_sum_tmpl(a, n); }

double reduce_sum(const double* a, int n) { return reduce_sum_tmpl(a, n); }

unsigned long long reduce_sum(const unsigned long long* a, int n) {
  return reduce_sum_tmpl(a, n);
}

//====================================================================//

template <class T>
void reduce_sum2_tmpl(T* __restrict__ h_ans, int hn, const T* __restrict__ v,
                      int nelem, int neach) {
  assert(hn <= neach);
  #pragma acc parallel loop gang deviceptr(v)
  for (int iv = 0; iv < hn; ++iv) {
    T ans = 0;
    #pragma acc loop vector reduction(+:ans)
    for (int ig = 0; ig < nelem; ++ig)
      ans += v[iv + ig * neach];
    h_ans[iv] = ans;
  }
}

void reduce_sum2(int* h_ans, int hn, const int* v, int nelem, int neach) {
  reduce_sum2_tmpl(h_ans, hn, v, nelem, neach);
}

void reduce_sum2(float* h_ans, int hn, const float* v, int nelem, int neach) {
  reduce_sum2_tmpl(h_ans, hn, v, nelem, neach);
}

void reduce_sum2(double* h_ans, int hn, const double* v, int nelem, int neach) {
  reduce_sum2_tmpl(h_ans, hn, v, nelem, neach);
}

void reduce_sum2(unsigned long long* h_ans, int hn, const unsigned long long* v,
                 int nelem, int neach) {
  reduce_sum2_tmpl(h_ans, hn, v, nelem, neach);
}

//====================================================================//

template <class T>
T dotprod_tmpl(const T* __restrict__ gpu_a, const T* __restrict__ gpu_b,
               int cpu_n) {
  T val = 0;
  #pragma acc parallel loop deviceptr(gpu_a,gpu_b) reduction(+:val)
  for (int i = 0; i < cpu_n; ++i)
    val += gpu_a[i] * gpu_b[i];
  return val;
}

float dotprod(const float* a, const float* b, int n) {
  return dotprod_tmpl(a, b, n);
}

double dotprod(const double* a, const double* b, int n) {
  return dotprod_tmpl(a, b, n);
}

//====================================================================//

template <class T>
void scale_array_tmpl(T* gpu_dst, T scal, int nelem) {
  #pragma acc parallel loop deviceptr(gpu_dst)
  for (int i = 0; i < nelem; ++i) {
    gpu_dst[i] *= scal;
  }
}

void scale_array(float* gpu_dst, float scal, int nelem) {
  scale_array_tmpl(gpu_dst, scal, nelem);
}

void scale_array(double* gpu_dst, double scal, int nelem) {
  scale_array_tmpl(gpu_dst, scal, nelem);
}
TINKER_NAMESPACE_END
