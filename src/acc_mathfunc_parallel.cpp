#include "acclib.h"
#include "deduce_ptr.h"
#include "mathfunc_parallel_acc.h"
#include "wait_queue.h"


TINKER_NAMESPACE_BEGIN
namespace platform {
namespace acc {
template <class T>
T reduce_sum(const T* gpu_a, size_t cpu_n, int sync)
{
   T val = 0;
   if (sync) {
      #pragma acc parallel loop independent\
                  deviceptr(gpu_a) reduction(+:val)
      for (size_t i = 0; i < cpu_n; ++i)
         val += gpu_a[i];
   } else {
      #pragma acc parallel loop independent async\
                  deviceptr(gpu_a) reduction(+:val)
      for (size_t i = 0; i < cpu_n; ++i)
         val += gpu_a[i];
   }
   return val;
}
template int reduce_sum(const int*, size_t, int);
template float reduce_sum(const float*, size_t, int);
template double reduce_sum(const double*, size_t, int);
template unsigned long long reduce_sum(const unsigned long long*, size_t, int);


template <class HT, size_t HN, class DPTR>
void reduce_sum2(HT (&restrict h_ans)[HN], DPTR restrict v, size_t nelem,
                 int sync)
{
   typedef typename deduce_ptr<DPTR>::type CONST_DT;
   typedef typename std::remove_const<CONST_DT>::type DT;
   static_assert(std::is_same<HT, DT>::value, "");

   constexpr size_t neach = deduce_ptr<DPTR>::n;
   static_assert(HN <= neach, "");

   for (size_t iv = 0; iv < HN; ++iv) {
      HT ans = 0;
      if (sync) {
         #pragma acc parallel loop independent\
                     deviceptr(v) reduction(+:ans)
         for (size_t ig = 0; ig < nelem; ++ig)
            ans += v[ig][iv];
      } else {
         #pragma acc parallel loop independent async\
                     deviceptr(v) reduction(+:ans)
         for (size_t ig = 0; ig < nelem; ++ig)
            ans += v[ig][iv];
      }
      h_ans[iv] = ans;
   }
}
template void reduce_sum2(float (&)[6], float (*)[8], size_t, int);
template void reduce_sum2(double (&)[6], double (*)[8], size_t, int);
template void reduce_sum2(unsigned long long (&)[6], unsigned long long (*)[8],
                          size_t, int);


template <class T>
T reduce_logic_or(const T* gpu_a, size_t cpu_n, int sync)
{
   T val = false;
   if (sync) {
      #pragma acc parallel loop independent\
                  deviceptr(gpu_a) reduction(||:val)
      for (size_t i = 0; i < cpu_n; ++i)
         val = (val || gpu_a[i]);
   } else {
      #pragma acc parallel loop independent async\
                  deviceptr(gpu_a) reduction(||:val)
      for (size_t i = 0; i < cpu_n; ++i)
         val = (val || gpu_a[i]);
   }
   return val;
}
template int reduce_logic_or(const int*, size_t, int);


template <class T>
T dotprod(const T* restrict gpu_a, const T* restrict gpu_b, size_t cpu_n)
{
   T val = 0;
   #pragma acc parallel loop independent deviceptr(gpu_a,gpu_b) reduction(+:val)
   for (size_t i = 0; i < cpu_n; ++i)
      val += gpu_a[i] * gpu_b[i];
   return val;
}
template float dotprod(const float*, const float*, size_t);
template double dotprod(const double*, const double*, size_t);


template <class T>
void dotprod(T* ans, const T* a, const T* b, int nelem, int sync)
{
   T v = 0;
   if (sync) {
      #pragma acc parallel loop independent deviceptr(a,b) reduction(+:v)
      for (size_t i = 0; i < nelem; ++i) {
         v += a[i] * b[i];
      }
      #pragma acc serial deviceptr(ans)
      {
         *ans = v;
      }
   } else {
      #pragma acc parallel loop independent async deviceptr(a,b) reduction(+:v)
      for (size_t i = 0; i < nelem; ++i) {
         v += a[i] * b[i];
      }
      #pragma acc serial async deviceptr(ans)
      {
         *ans = v;
      }
   }
}
template void dotprod(float*, const float*, const float*, int, int);
template void dotprod(double*, const double*, const double*, int, int);


template <class T>
void scale_array(T* gpu_dst, T scal, size_t nelem, DMFlag flag)
{
   if (flag & DMFlag::DEFAULT_Q) {
      #pragma acc parallel loop independent deviceptr(gpu_dst)
      for (size_t i = 0; i < nelem; ++i) {
         gpu_dst[i] *= scal;
      }
   } else {
      #pragma acc parallel loop independent async deviceptr(gpu_dst)
      for (size_t i = 0; i < nelem; ++i) {
         gpu_dst[i] *= scal;
      }
   }
   // if (flag & DMFlag::WAIT) {
   wait_queue(flag);
   // }
}
template void scale_array(float*, float, size_t, DMFlag);
template void scale_array(double*, double, size_t, DMFlag);
}
}
TINKER_NAMESPACE_END
