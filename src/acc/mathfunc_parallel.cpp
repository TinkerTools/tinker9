#include "glob.accasync.h"
#include "mathfunc_parallel_acc.h"
#include "tool/deduce_ptr.h"
#include <cassert>

namespace tinker {
template <class T>
T reduce_sum_acc(const T* gpu_a, size_t cpu_n, int queue)
{
   T val = 0;
   #pragma acc parallel loop independent async(queue)\
               deviceptr(gpu_a) copy(val) reduction(+:val)
   for (size_t i = 0; i < cpu_n; ++i)
      val += gpu_a[i];
   #pragma acc wait(queue)
   return val;
}
template int reduce_sum_acc(const int*, size_t, int);
template float reduce_sum_acc(const float*, size_t, int);
template double reduce_sum_acc(const double*, size_t, int);
template unsigned long long reduce_sum_acc(const unsigned long long*, size_t, int);

template <class HT, size_t HN, class DPTR>
void reduce_sum2_acc(HT (&restrict h_ans)[HN], DPTR restrict v, size_t nelem, int queue)
{
   typedef typename DeducePtr<DPTR>::type CONST_DT;
   typedef typename std::remove_const<CONST_DT>::type DT;
   static_assert(std::is_same<HT, DT>::value, "");

   constexpr size_t neach = DeducePtr<DPTR>::n;
   static_assert(HN <= neach, "");

   for (size_t iv = 0; iv < HN; ++iv) {
      HT ans = 0;
      #pragma acc parallel loop independent async(queue)\
                  deviceptr(v) copy(ans) reduction(+:ans)
      for (size_t ig = 0; ig < nelem; ++ig)
         ans += v[ig][iv];
      #pragma acc wait(queue)
      h_ans[iv] = ans;
   }
}
template void reduce_sum2_acc(float (&)[6], float (*)[8], size_t, int);
template void reduce_sum2_acc(double (&)[6], double (*)[8], size_t, int);
template void reduce_sum2_acc(unsigned long long (&)[6], unsigned long long (*)[8], size_t, int);

template <class T>
void reduce_sum_on_device_acc(T* dp_ans, const T* a, size_t nelem, int queue)
{
   static T ans1;
   static bool first = true;
   #pragma acc enter data if(first) async(queue) create(ans1)
   if (first)
      first = false;

   ans1 = 0;
   #pragma acc update async(queue) device(ans1)
   #pragma acc parallel loop independent async(queue) deviceptr(a)\
               reduction(+:ans1) present(ans1)
   for (size_t i = 0; i < nelem; ++i) {
      ans1 += a[i];
   }
   #pragma acc serial async(queue) present(ans1) deviceptr(dp_ans)
   {
      *dp_ans = ans1;
   }
}
template void reduce_sum_on_device_acc(int*, const int*, size_t, int);
template void reduce_sum_on_device_acc(float*, const float*, size_t, int);
template void reduce_sum_on_device_acc(double*, const double*, size_t, int);
template void reduce_sum_on_device_acc(unsigned long long*, const unsigned long long*, size_t, int);

template <class T, size_t HN, class DPTR>
void reduce_sum2_on_device_acc(T (&dref)[HN], DPTR v, size_t nelem, int queue)
{
   static T ans1;
   static bool first = true;
   #pragma acc enter data if(first) async(queue) create(ans1)
   if (first)
      first = false;

   T* dptr = &dref[0];
   for (size_t iv = 0; iv < HN; ++iv) {
      ans1 = 0;
      #pragma acc update async(queue) device(ans1)
      #pragma acc parallel loop independent async(queue) deviceptr(v)\
                  reduction(+:ans1) present(ans1)
      for (size_t ig = 0; ig < nelem; ++ig) {
         ans1 += v[ig][iv];
      }
      #pragma acc serial async(queue) deviceptr(dptr) present(ans1)
      {
         dptr[iv] = ans1;
      }
   }
}
template void reduce_sum2_on_device_acc(float (&)[6], float (*)[8], size_t, int);
template void reduce_sum2_on_device_acc(double (&)[6], double (*)[8], size_t, int);
template void reduce_sum2_on_device_acc(
   unsigned long long (&)[6], unsigned long long (*)[8], size_t, int);

template <class T>
T dotprod_acc(const T* restrict gpu_a, const T* restrict gpu_b, size_t cpu_n, int queue)
{
   T val = 0;
   #pragma acc parallel loop independent async(queue)\
               deviceptr(gpu_a,gpu_b) copy(val) reduction(+:val)
   for (size_t i = 0; i < cpu_n; ++i) {
      val += gpu_a[i] * gpu_b[i];
   }
   #pragma acc wait(queue)
   return val;
}
template float dotprod_acc(const float*, const float*, size_t, int);
template double dotprod_acc(const double*, const double*, size_t, int);

template <class T>
void dotprod_acc(T* ans, const T* a, const T* b, size_t nelem, int queue)
{
   static T ans1;
   static bool first = true;
   #pragma acc enter data if(first) async(queue) create(ans1)
   if (first)
      first = false;

   ans1 = 0;
   #pragma acc update async(queue) device(ans1)
   #pragma acc parallel loop independent async(queue) deviceptr(a,b)\
               reduction(+:ans1) present(ans1)
   for (size_t i = 0; i < nelem; ++i) {
      ans1 += a[i] * b[i];
   }
   #pragma acc serial async(queue) present(ans1) deviceptr(ans)
   {
      *ans = ans1;
   }
}
template void dotprod_acc(float*, const float*, const float*, size_t, int);
template void dotprod_acc(double*, const double*, const double*, size_t, int);

template <class T>
void scale_array_acc(T* gpu_dst, T scal, size_t nelem, int queue)
{
   #pragma acc parallel loop independent async(queue) deviceptr(gpu_dst)
   for (size_t i = 0; i < nelem; ++i) {
      gpu_dst[i] *= scal;
   }
}
template void scale_array_acc(float*, float, size_t, int);
template void scale_array_acc(double*, double, size_t, int);
}
