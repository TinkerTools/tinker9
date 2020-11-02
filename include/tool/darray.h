#pragma once
#include "glob.accasync.h"
#include "mathfunc.h"
#include "tool/deduce_ptr.h"
#include "tool/lpflag.h"
#include <vector>


namespace tinker {
void wait_for(int queue);
/**
 * \ingroup rc
 * Similar to OpenACC copyin, copies data from host to device.
 * \param dst     Device pointer.
 * \param src     Host pointer.
 * \param nbytes  Number of bytes.
 * \param flag    Kernel policy.
 * \see LPFlag
 */
void device_memory_copyin_bytes(void* dst, const void* src, size_t nbytes,
                                LPFlag flag);
/**
 * \ingroup rc
 * Similar to OpenACC copyout, copies data from device to host.
 * \param dst     Host pointer.
 * \param src     Device pointer.
 * \param nbytes  Number of bytes.
 * \param queue   OpenACC queue.
 */
void device_memory_copyout_bytes_async(void* dst, const void* src,
                                       size_t nbytes, int queue);
/**
 * \ingroup rc
 * Copies data between two pointers.
 * \note Different from OpenACC copy.
 * \param dst     Destination device pointer.
 * \param src     Source device pointer.
 * \param nbytes  Number of bytes.
 * \param flag    Kernel policy.
 * \see LPFlag
 */
void device_memory_copy_bytes(void* dst, const void* src, size_t nbytes,
                              LPFlag flag);
/**
 * \ingroup rc
 * Writes zero bytes on device.
 * \param dst     Device pointer.
 * \param nbytes  Number of bytes.
 * \param queue   OpenACC queue.
 */
void device_memory_zero_bytes_async(void* dst, size_t nbytes, int queue);
/**
 * \ingroup rc
 * Deallocates device pointer.
 * \param ptr  Device pointer.
 */
void device_memory_deallocate_bytes(void* ptr);
/**
 * \ingroup rc
 * Allocates device pointer.
 * \param pptr    Pointer to the device pointer.
 * \param nbytes  Number of bytes.
 */
void device_memory_allocate_bytes(void** pptr, size_t nbytes);
}


namespace tinker {
template <class T>
void device_memory_check_type()
{
   static_assert(std::is_enum<T>::value || std::is_integral<T>::value ||
                    std::is_floating_point<T>::value ||
                    std::is_trivial<T>::value,
                 "");
}


/**
 * \ingroup rc
 * Copies data to 1D array from host to device.
 * \param dst    Destination address.
 * \param src    Source address.
 * \param nelem  Number of elements to copy to the 1D device array.
 * \param flag   Kernel policy.
 * \see LPFlag
 */
template <class DT, class ST>
void device_memory_copyin_1d_array(DT* dst, const ST* src, size_t nelem,
                                   LPFlag flag)
{
   device_memory_check_type<DT>();
   device_memory_check_type<ST>();
   constexpr size_t ds = sizeof(DT); // device type
   constexpr size_t ss = sizeof(ST); // host type

   size_t size = ds * nelem;
   if (ds == ss) {
      device_memory_copyin_bytes(dst, src, size, flag);
   } else {
      std::vector<DT> buf(nelem);
      for (size_t i = 0; i < nelem; ++i)
         buf[i] = src[i];
      device_memory_copyin_bytes(dst, buf.data(), size, flag);
   }
}


/**
 * \ingroup rc
 * Copies data to 1D array from device to host.
 * \param dst    Destination address.
 * \param src    Source address.
 * \param nelem  Number of elements to copy to the 1D host array.
 * \param q      OpenACC queue.
 * \see LPFlag
 */
template <class DT, class ST>
void device_memory_copyout_1d_array(DT* dst, const ST* src, size_t nelem, int q)
{
   device_memory_check_type<DT>();
   device_memory_check_type<ST>();
   constexpr size_t ds = sizeof(DT); // host type
   constexpr size_t ss = sizeof(ST); // device type

   size_t size = ss * nelem;
   if (ds == ss) {
      device_memory_copyout_bytes_async(dst, src, size, q);
   } else {
      std::vector<ST> buf(nelem);
      device_memory_copyout_bytes_async(buf.data(), src, size, q);
      for (size_t i = 0; i < nelem; ++i)
         dst[i] = buf[i];
   }
}
}


namespace tinker {
/**
 * \ingroup rc
 * Device array.
 */
struct darray
{
   template <class T, size_t N>
   struct pointer;


   template <class T>
   struct pointer<T, 1>
   {
      typedef T* type;
   };


   template <class T, size_t N>
   struct pointer
   {
      static_assert(N > 1, "");
      typedef T (*type)[N];
   };


   template <class PTR>
   static typename deduce_ptr<PTR>::type* flatten(PTR p)
   {
      typedef typename deduce_ptr<PTR>::type T;
      return reinterpret_cast<T*>(p);
   }


   template <class PTR>
   static void allocate(size_t nelem, PTR* pp)
   {
      typedef typename deduce_ptr<PTR>::type T;
      constexpr size_t N = deduce_ptr<PTR>::n;
      device_memory_allocate_bytes(reinterpret_cast<void**>(pp),
                                   sizeof(T) * nelem * N);
   }


   template <class PTR, class... PTRS>
   static void allocate(size_t nelem, PTR* pp, PTRS... pps)
   {
      allocate(nelem, pp);
      allocate(nelem, pps...);
   }


   template <class PTR>
   static void deallocate(PTR p)
   {
      device_memory_deallocate_bytes(flatten(p));
   }


   template <class PTR, class... PTRS>
   static void deallocate(PTR p, PTRS... ps)
   {
      deallocate(p);
      deallocate(ps...);
   }


   template <class PTR>
   static void zero(int q, size_t nelem, PTR p)
   {
      typedef typename deduce_ptr<PTR>::type T;
      constexpr size_t N = deduce_ptr<PTR>::n;
      device_memory_zero_bytes_async(flatten(p), sizeof(T) * nelem * N, q);
   }


   template <class PTR, class... PTRS>
   static void zero(int q, size_t nelem, PTR p, PTRS... ps)
   {
      zero(q, nelem, p);
      zero(q, nelem, ps...);
   }


   template <class PTR, class U>
   static void copyin(LPFlag flag, size_t nelem, PTR dst, const U* src)
   {
      constexpr size_t N = deduce_ptr<PTR>::n;
      device_memory_copyin_1d_array(flatten(dst), flatten(src), nelem * N,
                                    flag);
   }


   template <class U, class PTR>
   static void copyout(size_t nelem, U* dst, const PTR src, int q)
   {
      constexpr size_t N = deduce_ptr<PTR>::n;
      device_memory_copyout_1d_array(flatten(dst), flatten(src), nelem * N, q);
   }


   template <class PTR, class U>
   static void copy(LPFlag flag, size_t nelem, PTR dst, const U* src)
   {
      constexpr size_t N = deduce_ptr<PTR>::n;
      using DT = typename deduce_ptr<PTR>::type;
      using ST = typename deduce_ptr<U*>::type;
      static_assert(std::is_same<DT, ST>::value, "");
      size_t size = N * sizeof(ST) * nelem;
      device_memory_copy_bytes(flatten(dst), flatten(src), size, flag);
   }


   template <class DT, class ST>
   static void copyin2(LPFlag flag, size_t idx0, size_t ndim, size_t nelem,
                       DT dst, const ST src)
   {
      static_assert(deduce_ptr<DT>::n == 1, "");
      static_assert(deduce_ptr<ST>::n == 1, "");
      typedef typename deduce_ptr<DT>::type T;
      std::vector<T> buf(nelem);
      for (size_t i = 0; i < nelem; ++i)
         buf[i] = src[ndim * i + idx0];
      copyin(flag, nelem, dst, buf.data());
   }


   template <class DT, class ST>
   static void copyout2(LPFlag flag, size_t idx0, size_t ndim, size_t nelem,
                        DT dst, const ST src)
   {
      static_assert(deduce_ptr<DT>::n == 1, "");
      static_assert(deduce_ptr<ST>::n == 1, "");
      typedef typename deduce_ptr<ST>::type T;
      std::vector<T> buf(nelem);
      copyout(flag, nelem, buf.data(), src);
      for (size_t i = 0; i < nelem; ++i)
         dst[ndim * i + idx0] = buf[i];
   }


   template <class PTR, class PTR2>
   static typename deduce_ptr<PTR>::type dot(LPFlag flag, size_t nelem,
                                             const PTR ptr, const PTR2 b)
   {
      typedef typename deduce_ptr<PTR>::type T;
      constexpr size_t N = deduce_ptr<PTR>::n;
      typedef typename deduce_ptr<PTR2>::type T2;
      static_assert(std::is_same<T, T2>::value, "");
      return parallel::dotprod(flatten(ptr), flatten(b), nelem * N, flag);
   }


   template <class ANS, class PTR, class PTR2>
   static void dot(LPFlag flag, int nelem, ANS ans, const PTR ptr,
                   const PTR2 ptr2)
   {
      typedef typename deduce_ptr<PTR>::type T;
      constexpr size_t N = deduce_ptr<PTR>::n;
      typedef typename deduce_ptr<PTR2>::type T2;
      static_assert(std::is_same<T, T2>::value, "");
      typedef typename deduce_ptr<ANS>::type TA;
      static_assert(std::is_same<T, TA>::value, "");
      parallel::dotprod(ans, flatten(ptr), flatten(ptr2), nelem * N, flag);
   }


   template <class FLT, class PTR>
   static void scale(LPFlag flag, size_t nelem, FLT scal, PTR ptr)
   {
      constexpr size_t N = deduce_ptr<PTR>::n;
      parallel::scale_array(flatten(ptr), scal, nelem * N, flag);
   }


   template <class FLT, class PTR, class... PTRS>
   static void scale(LPFlag flag, size_t nelem, FLT scal, PTR ptr, PTRS... ptrs)
   {
      scale(flag, nelem, scal, ptr);
      scale(flag, nelem, scal, ptrs...);
   }
};


/**
 * \ingroup rc
 * Based on the template parameters `T` and `N`, this type is either
 * defined to `T*` or `T(*)[N]` when `N` is greater than 1.
 * `N` is set to 1 by default.
 */
template <class T, size_t N = 1>
using pointer = typename darray::pointer<T, N>::type;
}
