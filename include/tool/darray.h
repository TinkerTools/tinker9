#pragma once
#include "math/const.h"
#include "math/libfunc.h"
#include "math/parallel.h"
#include "tool/accasync.h"
#include "tool/ptrtrait.h"

#include <vector>

namespace tinker {
/// \ingroup rc
/// \brief Similar to OpenACC wait and CUDA stream synchronize.
void waitFor(int queue ///< OpenACC queue.
);

/// \ingroup rc
/// \brief Similar to OpenACC async copyin, copies data from host to device.
void deviceMemoryCopyinBytesAsync(void* dst,       ///< Device pointer.
                                  const void* src, ///< Host pointer.
                                  size_t nbytes,   ///< Number of bytes.
                                  int queue        ///< OpenACC queue.
);

/// \ingroup rc
/// \brief Similar to OpenACC async copyout, copies data from device to host.
void deviceMemoryCopyoutBytesAsync(void* dst,       ///< Host pointer.
                                   const void* src, ///< Device pointer.
                                   size_t nbytes,   ///< Number of bytes.
                                   int queue        ///< OpenACC queue.
);

/// \ingroup rc
/// \brief Copies data between two pointers on device.
/// \note Different from OpenACC copy.
void deviceMemoryCopyBytesAsync(void* dst,       ///< Destination device pointer.
                                const void* src, ///< Source device pointer.
                                size_t nbytes,   ///< Number of bytes.
                                int queue        ///< OpenACC queue.
);

/// \ingroup rc
/// \brief Writes zero bytes on device.
void deviceMemoryZeroBytesAsync(void* dst,     ///< Device pointer.
                                size_t nbytes, ///< Number of bytes.
                                int queue      ///< OpenACC queue.
);

/// \ingroup rc
/// \brief Deallocates device pointer.
void deviceMemoryDeallocate(void* ptr ///< Device pointer.
);

/// \ingroup rc
/// \brief Allocates device pointer.
void deviceMemoryAllocateBytes(void** pptr,  ///< Pointer to the device pointer.
                               size_t nbytes ///< Number of bytes.
);
}

namespace tinker {
inline namespace v1 {
/// \ingroup rc
/// \brief Sanity check.
template <class T>
void deviceMemoryCheckType()
{
   static_assert(std::is_enum<T>::value || std::is_integral<T>::value || std::is_floating_point<T>::value || std::is_trivial<T>::value, "");
}

/// \ingroup rc
/// \brief Copies data to 1D array, host to device.
template <class DT, class ST>
void deviceMemoryCopyin1dArray(DT* dst,       ///< Destination address.
                               const ST* src, ///< Source address.
                               size_t nelem,  ///< Number of elements to copy to the 1D device array.
                               int q          ///< OpenACC queue.
)
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
      waitFor(q);
   }
}

/// \ingroup rc
/// \brief Copies data to 1D array, device to host.
template <class DT, class ST>
void deviceMemoryCopyout1dArray(DT* dst,       ///< Destination address.
                                const ST* src, ///< Source address.
                                size_t nelem,  ///< Number of elements to copy to the 1D host array.
                                int q          ///< OpenACC queue.
)
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
      waitFor(q);
      for (size_t i = 0; i < nelem; ++i)
         dst[i] = buf[i];
   }
}
}
}

namespace tinker {
/// \ingroup rc
/// \brief Device array.
class darray
{
private:
   template <class PTR>
   static typename PtrTrait<PTR>::type* flatten(PTR p)
   {
      typedef typename PtrTrait<PTR>::type T;
      return reinterpret_cast<T*>(p);
   }

public:
   template <class PTR>
   static void allocate(size_t nelem, PTR* pp)
   {
      typedef typename PtrTrait<PTR>::type T;
      constexpr size_t N = PtrTrait<PTR>::n;
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
      typedef typename PtrTrait<PTR>::type T;
      constexpr size_t N = PtrTrait<PTR>::n;
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
      constexpr size_t N = PtrTrait<PTR>::n;
      deviceMemoryCopyin1dArray(flatten(dst), flatten(src), nelem * N, q);
   }

   template <class U, class PTR>
   static void copyout(int q, size_t nelem, U* dst, const PTR src)
   {
      constexpr size_t N = PtrTrait<PTR>::n;
      deviceMemoryCopyout1dArray(flatten(dst), flatten(src), nelem * N, q);
   }

   /// \brief Copies data across two device memory pointers.
   template <class PTR, class U>
   static void copy(int q, size_t nelem, PTR dst, const U* src)
   {
      constexpr size_t N = PtrTrait<PTR>::n;
      using DT = typename PtrTrait<PTR>::type;
      using ST = typename PtrTrait<U*>::type;
      static_assert(std::is_same<DT, ST>::value, "");
      size_t size = N * sizeof(ST) * nelem;
      deviceMemoryCopyBytesAsync(flatten(dst), flatten(src), size, q);
   }

   /// \brief Calculates the dot product and returns the answer to the host.
   template <class PTR, class PTR2>
   static typename PtrTrait<PTR>::type dotThenReturn(int q, size_t nelem, const PTR ptr, const PTR2 b)
   {
      typedef typename PtrTrait<PTR>::type T;
      constexpr size_t N = PtrTrait<PTR>::n;
      typedef typename PtrTrait<PTR2>::type T2;
      static_assert(std::is_same<T, T2>::value, "");
      return dotProd(flatten(ptr), flatten(b), nelem * N, q);
   }

   /// \brief Calculates the dot product and saves the answer to pointer `ans`.
   template <class ANS, class PTR, class PTR2>
   static void dot(int q, size_t nelem, ANS ans, const PTR ptr, const PTR2 ptr2)
   {
      typedef typename PtrTrait<PTR>::type T;
      constexpr size_t N = PtrTrait<PTR>::n;
      typedef typename PtrTrait<PTR2>::type T2;
      static_assert(std::is_same<T, T2>::value, "");
      typedef typename PtrTrait<ANS>::type TA;
      static_assert(std::is_same<T, TA>::value, "");
      dotProd(ans, flatten(ptr), flatten(ptr2), nelem * N, q);
   }

   template <class FLT, class PTR>
   static void scale(int q, size_t nelem, FLT scal, PTR ptr)
   {
      constexpr size_t N = PtrTrait<PTR>::n;
      scaleArray(flatten(ptr), scal, nelem * N, q);
   }

   template <class FLT, class PTR, class... PTRS>
   static void scale(int q, size_t nelem, FLT scal, PTR ptr, PTRS... ptrs)
   {
      scale(q, nelem, scal, ptr);
      scale(q, nelem, scal, ptrs...);
   }
};
}
