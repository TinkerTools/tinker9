#pragma once
#include "glob.accasync.h"
#include "mathfunc.h"
#include "tool/deduce_ptr.h"
#include <vector>

namespace tinker {
/**
 * \ingroup rc
 * Similar to OpenACC wait and CUDA stream synchronize.
 * \param queue  OpenACC queue.
 */
void wait_for(int queue);
/**
 * \ingroup rc
 * Similar to OpenACC async copyin, copies data from host to device.
 * \param dst     Device pointer.
 * \param src     Host pointer.
 * \param nbytes  Number of bytes.
 * \param queue   OpenACC queue.
 */
void deviceMemoryCopyinBytesAsync(void* dst, const void* src, size_t nbytes, int queue);
/**
 * \ingroup rc
 * Similar to OpenACC async copyout, copies data from device to host.
 * \param dst     Host pointer.
 * \param src     Device pointer.
 * \param nbytes  Number of bytes.
 * \param queue   OpenACC queue.
 */
void deviceMemoryCopyoutBytesAsync(void* dst, const void* src, size_t nbytes, int queue);
/**
 * \ingroup rc
 * Copies data between two pointers.
 * \note Different from OpenACC copy.
 * \param dst     Destination device pointer.
 * \param src     Source device pointer.
 * \param nbytes  Number of bytes.
 * \param queue   OpenACC queue.
 */
void deviceMemoryCopyBytesAsync(void* dst, const void* src, size_t nbytes, int queue);
/**
 * \ingroup rc
 * Writes zero bytes on device.
 * \param dst     Device pointer.
 * \param nbytes  Number of bytes.
 * \param queue   OpenACC queue.
 */
void deviceMemoryZeroBytesAsync(void* dst, size_t nbytes, int queue);
/**
 * \ingroup rc
 * Deallocates device pointer.
 * \param ptr  Device pointer.
 */
void deviceMemoryDeallocate(void* ptr);
/**
 * \ingroup rc
 * Allocates device pointer.
 * \param pptr    Pointer to the device pointer.
 * \param nbytes  Number of bytes.
 */
void deviceMemoryAllocateBytes(void** pptr, size_t nbytes);
}

namespace tinker {
template <class T>
void deviceMemoryCheckType()
{
   static_assert(std::is_enum<T>::value || std::is_integral<T>::value ||
         std::is_floating_point<T>::value || std::is_trivial<T>::value,
      "");
}

/**
 * \ingroup rc
 * Copies data to 1D array, host to device.
 * \param dst    Destination address.
 * \param src    Source address.
 * \param nelem  Number of elements to copy to the 1D device array.
 * \param q      OpenACC queue.
 */
template <class DT, class ST>
void deviceMemoryCopyin1dArray(DT* dst, const ST* src, size_t nelem, int q)
{
   deviceMemoryCheckType<DT>();
   deviceMemoryCheckType<ST>();
   constexpr size_t ds = sizeof(DT); // device type
   constexpr size_t ss = sizeof(ST); // host type

   size_t size = ds * nelem;
   if (ds == ss) {
      deviceMemoryCopyinBytesAsync(dst, src, size, q);
   } else {
      std::vector<DT> buf(nelem);
      for (size_t i = 0; i < nelem; ++i)
         buf[i] = src[i];
      deviceMemoryCopyinBytesAsync(dst, buf.data(), size, q);
      wait_for(q);
   }
}

/**
 * \ingroup rc
 * Copies data to 1D array, device to host.
 * \param dst    Destination address.
 * \param src    Source address.
 * \param nelem  Number of elements to copy to the 1D host array.
 * \param q      OpenACC queue.
 * \see LPFlag
 */
template <class DT, class ST>
void deviceMemoryCopyout1dArray(DT* dst, const ST* src, size_t nelem, int q)
{
   deviceMemoryCheckType<DT>();
   deviceMemoryCheckType<ST>();
   constexpr size_t ds = sizeof(DT); // host type
   constexpr size_t ss = sizeof(ST); // device type

   size_t size = ss * nelem;
   if (ds == ss) {
      deviceMemoryCopyoutBytesAsync(dst, src, size, q);
   } else {
      std::vector<ST> buf(nelem);
      deviceMemoryCopyoutBytesAsync(buf.data(), src, size, q);
      wait_for(q);
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
   static typename DeducePtr<PTR>::type* flatten(PTR p)
   {
      typedef typename DeducePtr<PTR>::type T;
      return reinterpret_cast<T*>(p);
   }

   template <class PTR>
   static void allocate(size_t nelem, PTR* pp)
   {
      typedef typename DeducePtr<PTR>::type T;
      constexpr size_t N = DeducePtr<PTR>::n;
      deviceMemoryAllocateBytes(reinterpret_cast<void**>(pp), sizeof(T) * nelem * N);
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
      deviceMemoryDeallocate(flatten(p));
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
      typedef typename DeducePtr<PTR>::type T;
      constexpr size_t N = DeducePtr<PTR>::n;
      deviceMemoryZeroBytesAsync(flatten(p), sizeof(T) * nelem * N, q);
   }

   template <class PTR, class... PTRS>
   static void zero(int q, size_t nelem, PTR p, PTRS... ps)
   {
      zero(q, nelem, p);
      zero(q, nelem, ps...);
   }

   template <class PTR, class U>
   static void copyin(int q, size_t nelem, PTR dst, const U* src)
   {
      constexpr size_t N = DeducePtr<PTR>::n;
      deviceMemoryCopyin1dArray(flatten(dst), flatten(src), nelem * N, q);
   }

   template <class U, class PTR>
   static void copyout(int q, size_t nelem, U* dst, const PTR src)
   {
      constexpr size_t N = DeducePtr<PTR>::n;
      deviceMemoryCopyout1dArray(flatten(dst), flatten(src), nelem * N, q);
   }

   template <class PTR, class U>
   static void copy(int q, size_t nelem, PTR dst, const U* src)
   {
      constexpr size_t N = DeducePtr<PTR>::n;
      using DT = typename DeducePtr<PTR>::type;
      using ST = typename DeducePtr<U*>::type;
      static_assert(std::is_same<DT, ST>::value, "");
      size_t size = N * sizeof(ST) * nelem;
      deviceMemoryCopyBytesAsync(flatten(dst), flatten(src), size, q);
   }

   template <class PTR, class PTR2>
   static typename DeducePtr<PTR>::type dot_wait(int q, size_t nelem, const PTR ptr, const PTR2 b)
   {
      typedef typename DeducePtr<PTR>::type T;
      constexpr size_t N = DeducePtr<PTR>::n;
      typedef typename DeducePtr<PTR2>::type T2;
      static_assert(std::is_same<T, T2>::value, "");
      return dotprod(flatten(ptr), flatten(b), nelem * N, q);
   }

   template <class ANS, class PTR, class PTR2>
   static void dot(int q, size_t nelem, ANS ans, const PTR ptr, const PTR2 ptr2)
   {
      typedef typename DeducePtr<PTR>::type T;
      constexpr size_t N = DeducePtr<PTR>::n;
      typedef typename DeducePtr<PTR2>::type T2;
      static_assert(std::is_same<T, T2>::value, "");
      typedef typename DeducePtr<ANS>::type TA;
      static_assert(std::is_same<T, TA>::value, "");
      dotprod(ans, flatten(ptr), flatten(ptr2), nelem * N, q);
   }

   template <class FLT, class PTR>
   static void scale(int q, size_t nelem, FLT scal, PTR ptr)
   {
      constexpr size_t N = DeducePtr<PTR>::n;
      scale_array(flatten(ptr), scal, nelem * N, q);
   }

   template <class FLT, class PTR, class... PTRS>
   static void scale(int q, size_t nelem, FLT scal, PTR ptr, PTRS... ptrs)
   {
      scale(q, nelem, scal, ptr);
      scale(q, nelem, scal, ptrs...);
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
