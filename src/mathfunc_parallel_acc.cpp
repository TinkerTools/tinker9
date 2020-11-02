#include "mathfunc_parallel_acc.h"
#include "glob.accasync.h"
#include "tool/acclib.h"
#include "tool/deduce_ptr.h"
#include <cassert>


namespace tinker {
template <class T>
T reduce_sum_acc(const T* gpu_a, size_t cpu_n, LPFlag flag)
{
   T val = 0;
   if (flag & LPFlag::DEFAULT_Q) {
      #pragma acc parallel loop independent\
                  deviceptr(gpu_a) reduction(+:val)
      for (size_t i = 0; i < cpu_n; ++i)
         val += gpu_a[i];
   } else {
      // This OpenACC directive was
      // acc parallel loop independent async deviceptr(gpu_a) reduction(+:val)
      // but kept receiving segmentation fault. Luckily, this code is easy
      // enough for compiler to recognize it is reduction.
      #pragma acc parallel loop async deviceptr(gpu_a)
      for (size_t i = 0; i < cpu_n; ++i)
         val += gpu_a[i];
   }
   // implicit OpenACC wait
   assert(flag & LPFlag::WAIT);
   return val;
}
template int reduce_sum_acc(const int*, size_t, LPFlag);
template float reduce_sum_acc(const float*, size_t, LPFlag);
template double reduce_sum_acc(const double*, size_t, LPFlag);
template unsigned long long reduce_sum_acc(const unsigned long long*, size_t,
                                           LPFlag);


template <class HT, size_t HN, class DPTR>
void reduce_sum2_acc(HT (&restrict h_ans)[HN], DPTR restrict v, size_t nelem,
                     LPFlag flag)
{
   typedef typename deduce_ptr<DPTR>::type CONST_DT;
   typedef typename std::remove_const<CONST_DT>::type DT;
   static_assert(std::is_same<HT, DT>::value, "");

   constexpr size_t neach = deduce_ptr<DPTR>::n;
   static_assert(HN <= neach, "");

   bool sync = flag & LPFlag::DEFAULT_Q;
   for (size_t iv = 0; iv < HN; ++iv) {
      HT ans = 0;
      if (sync) {
         #pragma acc parallel loop independent\
                     deviceptr(v) reduction(+:ans)
         for (size_t ig = 0; ig < nelem; ++ig)
            ans += v[ig][iv];
      } else {
         // see reduce_sum()
         #pragma acc parallel loop async deviceptr(v)
         for (size_t ig = 0; ig < nelem; ++ig)
            ans += v[ig][iv];
      }
      // implicit OpenACC wait
      assert(flag & LPFlag::WAIT);
      h_ans[iv] = ans;
   }
}
template void reduce_sum2_acc(float (&)[6], float (*)[8], size_t, LPFlag);
template void reduce_sum2_acc(double (&)[6], double (*)[8], size_t, LPFlag);
template void reduce_sum2_acc(unsigned long long (&)[6],
                              unsigned long long (*)[8], size_t, LPFlag);


template <class T>
T reduce_logic_or_acc(const T* gpu_a, size_t cpu_n, LPFlag flag)
{
   T val = false;
   if (flag & LPFlag::DEFAULT_Q) {
      #pragma acc parallel loop independent\
                  deviceptr(gpu_a) reduction(||:val)
      for (size_t i = 0; i < cpu_n; ++i)
         val = (val || gpu_a[i]);
   } else {
      // see reduce_sum()
      #pragma acc parallel loop async deviceptr(gpu_a)
      for (size_t i = 0; i < cpu_n; ++i)
         val = (val || gpu_a[i]);
   }
   // implicit OpenACC wait
   assert(flag & LPFlag::WAIT);
   return val;
}
template int reduce_logic_or_acc(const int*, size_t, LPFlag);


template <class T>
T dotprod_acc(const T* restrict gpu_a, const T* restrict gpu_b, size_t cpu_n,
              LPFlag flag)
{
   T val = 0;
   if (flag & LPFlag::DEFAULT_Q) {
      #pragma acc parallel loop independent\
                  deviceptr(gpu_a,gpu_b) reduction(+:val)
      for (size_t i = 0; i < cpu_n; ++i)
         val += gpu_a[i] * gpu_b[i];
   } else {
      // see reduce_sum()
      #pragma acc parallel loop async deviceptr(gpu_a,gpu_b)
      for (size_t i = 0; i < cpu_n; ++i)
         val += gpu_a[i] * gpu_b[i];
   }
   // implicit OpenACC wait
   assert(flag & LPFlag::WAIT);
   return val;
}
template float dotprod_acc(const float*, const float*, size_t, LPFlag);
template double dotprod_acc(const double*, const double*, size_t, LPFlag);


template <class T>
void dotprod_acc(T* ans, const T* a, const T* b, int nelem, LPFlag flag)
{
   bool sync = flag & LPFlag::DEFAULT_Q;
   // T v = 0;
   if (sync) {
      #pragma acc serial deviceptr(ans)
      {
         *ans = 0;
      }
      #pragma acc parallel loop deviceptr(ans,a,b)
      for (int i = 0; i < nelem; ++i) {
         *ans += a[i] * b[i];
      }
   } else {
      #pragma acc serial async deviceptr(ans)
      {
         *ans = 0;
      }
      #pragma acc parallel loop async deviceptr(ans,a,b)
      for (int i = 0; i < nelem; ++i) {
         *ans += a[i] * b[i];
      }
   }
   // if (flag & LPFlag::WAIT) {
   wait_queue(flag);
   // }
}
template void dotprod_acc(float*, const float*, const float*, int, LPFlag);
template void dotprod_acc(double*, const double*, const double*, int, LPFlag);


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
