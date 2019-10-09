#pragma once

#include "deduce_ptr.h"

TINKER_NAMESPACE_BEGIN
/// \defgroup mem Memory Management
/// \ingroup gvar

void device_memory_copyin_bytes(void* dst, const void* src, size_t nbytes);
void device_memory_copyout_bytes(void* dst, const void* src, size_t nbytes);
void device_memory_copy_bytes(void* dst, const void* src, size_t nbytes);
void device_memory_zero_bytes(void* dst, size_t nbytes);
void device_memory_deallocate_bytes(void* ptr);
void device_memory_allocate_bytes(void** pptr, size_t nbytes);
TINKER_NAMESPACE_END

#include <stdexcept>
#include <type_traits>
#include <vector>

TINKER_NAMESPACE_BEGIN
struct DeviceMemory
{
   static void copyin_bytes(void* dst, const void* src, size_t nbytes)
   {
      device_memory_copyin_bytes(dst, src, nbytes);
   }

   static void copyout_bytes(void* dst, const void* src, size_t nbytes)
   {
      device_memory_copyout_bytes(dst, src, nbytes);
   }

   static void copy_bytes(void* dst, const void* src, size_t nbytes)
   {
      device_memory_copyout_bytes(dst, src, nbytes);
   }

   static void zero_bytes(void* dst, size_t nbytes)
   {
      device_memory_zero_bytes(dst, nbytes);
   }

   static void deallocate_bytes(void* ptr)
   {
      device_memory_deallocate_bytes(ptr);
   }

   static void allocate_bytes(void** pptr, size_t nbytes)
   {
      device_memory_allocate_bytes(pptr, nbytes);
      if (nbytes && *pptr == nullptr)
         throw std::bad_alloc();
   }

   //====================================================================//

   template <
      class T,
      class = typename std::enable_if<!std::is_same<T, void>::value>::type>
   static void allocate_bytes(T** pptr, size_t nbytes)
   {
      return allocate_bytes(reinterpret_cast<void**>(pptr), nbytes);
   }

   //====================================================================//

   struct CopyIn
   {
      void operator()(void* dst, const void* src, size_t nbytes)
      {
         copyin_bytes(dst, src, nbytes);
      }
   };

   struct CopyOut
   {
      void operator()(void* dst, const void* src, size_t nbytes)
      {
         copyout_bytes(dst, src, nbytes);
      }
   };

   struct Copy
   {
      void operator()(void* dst, const void* src, size_t nbytes)
      {
         copy_bytes(dst, src, nbytes);
      }
   };

   struct Zero
   {
      void operator()(void* ptr, size_t nbytes)
      {
         zero_bytes(ptr, nbytes);
      }
   };

   struct Deallocate
   {
      void operator()(void* ptr)
      {
         deallocate_bytes(ptr);
      }
   };

   struct Allocate
   {
      void operator()(void** pptr, size_t nbytes)
      {
         allocate_bytes(pptr, nbytes);
      }
   };

   //====================================================================//

   template <class T>
   static void check_type()
   {
      static_assert(std::is_enum<T>::value || std::is_integral<T>::value ||
                       std::is_floating_point<T>::value ||
                       std::is_trivial<T>::value,
                    "");
   }

   //====================================================================//

   template <class DT, class ST>
   static void copyin_array(DT* dst, const ST* src, size_t nelem)
   {
      check_type<DT>();
      check_type<ST>();
      constexpr size_t ds = sizeof(DT); // device type
      constexpr size_t ss = sizeof(ST); // host type

      size_t size = ds * nelem;
      if_constexpr(ds == ss)
      {
         copyin_bytes(dst, src, size);
      }
      else
      {
         std::vector<DT> buf(nelem);
         for (size_t i = 0; i < nelem; ++i)
            buf[i] = src[i];
         copyin_bytes(dst, buf.data(), size);
      }
   }

   template <class DT, class ST>
   static void copyout_array(DT* dst, const ST* src, size_t nelem)
   {
      check_type<DT>();
      check_type<ST>();
      constexpr size_t ds = sizeof(DT); // host type
      constexpr size_t ss = sizeof(ST); // device type

      size_t size = ss * nelem;
      if_constexpr(ds == ss)
      {
         copyout_bytes(dst, src, size);
      }
      else
      {
         std::vector<ST> buf(nelem);
         copyout_bytes(buf.data(), src, size);
         for (size_t i = 0; i < nelem; ++i)
            dst[i] = buf[i];
      }
   }

   template <class DT, class ST>
   static void copy_array(DT* dst, const ST* src, size_t nelem)
   {
      check_type<DT>();
      check_type<ST>();
      static_assert(std::is_same<DT, ST>::value, "");
      size_t size = sizeof(ST) * nelem;
      copy_bytes(dst, src, size);
   }

   template <class T>
   static void zero_array(T* dst, size_t nelem)
   {
      check_type<T>();
      size_t size = sizeof(T) * nelem;
      zero_bytes(dst, size);
   }
};
TINKER_NAMESPACE_END

#include "mathfunc.h"
TINKER_NAMESPACE_BEGIN
template <class Op>
struct DeviceArray
{
   template <class T, size_t N>
   struct ptr;

   template <class T>
   struct ptr<T, 1>
   {
      typedef T* type;
   };

   template <class T, size_t N>
   struct ptr
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

   //====================================================================//

   template <class PTR>
   static void allocate(size_t nelem, PTR* pp)
   {
      typedef typename deduce_ptr<PTR>::type T;
      constexpr size_t N = deduce_ptr<PTR>::n;
      Op a;
      a.allocate_bytes(reinterpret_cast<void**>(pp), sizeof(T) * nelem * N);
   }

   template <class PTR>
   static void deallocate(PTR p)
   {
      Op a;
      a.deallocate_bytes(flatten(p));
   }

   template <class PTR, class... PTRS>
   static void allocate(size_t nelem, PTR* pp, PTRS... pps)
   {
      allocate(nelem, pp);
      allocate(nelem, pps...);
   }

   template <class PTR, class... PTRS>
   static void deallocate(PTR p, PTRS... ps)
   {
      deallocate(p);
      deallocate(ps...);
   }

   //====================================================================//

   template <class PTR>
   static void zero(size_t nelem, PTR p)
   {
      typedef typename deduce_ptr<PTR>::type T;
      constexpr size_t N = deduce_ptr<PTR>::n;
      Op a;
      a.zero_bytes(flatten(p), sizeof(T) * nelem * N);
   }

   template <class PTR, class... PTRS>
   static void zero(size_t nelem, PTR p, PTRS... ps)
   {
      zero(nelem, p);
      zero(nelem, ps...);
   }

   //====================================================================//

   template <class PTR, class U>
   static void copyin(size_t nelem, PTR dst, const U* src)
   {
      typedef typename deduce_ptr<PTR>::type T;
      constexpr size_t N = deduce_ptr<PTR>::n;
      Op a;
      a.copyin_array(flatten(dst), flatten(src), nelem * N);
   }

   template <class U, class PTR>
   static void copyout(size_t nelem, U* dst, const PTR src)
   {
      typedef typename deduce_ptr<PTR>::type T;
      constexpr size_t N = deduce_ptr<PTR>::n;
      Op a;
      a.copyout_array(flatten(dst), flatten(src), nelem * N);
   }

   template <class PTR, class U>
   static void copy(size_t nelem, PTR dst, const U* src)
   {
      typedef typename deduce_ptr<PTR>::type T;
      constexpr size_t N = deduce_ptr<PTR>::n;
      Op a;
      a.copy_array(flatten(dst), flatten(src), nelem * N);
   }

   //====================================================================//

   template <class DT, class ST>
   static void copyin2(size_t idx0, size_t ndim, size_t nelem, DT dst,
                       const ST src)
   {
      static_assert(deduce_ptr<DT>::n == 1, "");
      static_assert(deduce_ptr<ST>::n == 1, "");
      typedef typename deduce_ptr<DT>::type T;
      std::vector<T> buf(nelem);
      for (size_t i = 0; i < nelem; ++i)
         buf[i] = src[ndim * i + idx0];
      copyin(nelem, dst, buf.data());
   }

   template <class DT, class ST>
   static void copyout2(size_t idx0, size_t ndim, size_t nelem, DT dst,
                        const ST src)
   {
      static_assert(deduce_ptr<DT>::n == 1, "");
      static_assert(deduce_ptr<ST>::n == 1, "");
      typedef typename deduce_ptr<ST>::type T;
      std::vector<T> buf(nelem);
      copyout(nelem, buf.data(), src);
      for (size_t i = 0; i < nelem; ++i)
         dst[ndim * i + idx0] = buf[i];
   }

   //====================================================================//

   template <class PTR, class PTR2>
   static typename deduce_ptr<PTR>::type dot(size_t nelem, const PTR ptr,
                                             const PTR2 b)
   {
      typedef typename deduce_ptr<PTR>::type T;
      constexpr size_t N = deduce_ptr<PTR>::n;
      typedef typename deduce_ptr<PTR2>::type T2;
      static_assert(std::is_same<T, T2>::value, "");
      return parallel::dotprod(flatten(ptr), flatten(b), nelem * N);
   }

   //====================================================================//

   template <class FLT, class PTR>
   static void scale(size_t nelem, FLT scal, PTR ptr)
   {
      typedef typename deduce_ptr<PTR>::type T;
      constexpr size_t N = deduce_ptr<PTR>::n;
      parallel::scale_array(flatten(ptr), scal, nelem * N);
   }

   template <class FLT, class PTR, class... PTRS>
   static void scale(size_t nelem, FLT scal, PTR ptr, PTRS... ptrs)
   {
      scale(nelem, scal, ptr);
      scale(nelem, scal, ptrs...);
   }
};
TINKER_NAMESPACE_END
