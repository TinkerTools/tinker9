#ifndef TINKER_DEV_ARRAY_H_
#define TINKER_DEV_ARRAY_H_

#include "dev_allocator.h"
#include "mathfunc.h"

TINKER_NAMESPACE_BEGIN
template <template <class> class Allocator>
class DeviceArray {
private:
  template <class PTR>
  struct deduce;

  template <class T>
  struct deduce<T*> {
    typedef T type;
    static constexpr size_t N = 1;
  };

  template <class T, size_t N1>
  struct deduce<T (*)[N1]> {
    typedef T type;
    static constexpr size_t N = N1;
  };

public:
  template <class T, size_t N = 1>
  struct ptr;

  template <class T>
  struct ptr<T> {
    typedef T* type;
  };

  template <class T, size_t N>
  struct ptr {
    static_assert(N > 1, "");
    typedef T (*type)[N];
  };

  template <class PTR>
  static typename deduce<PTR>::type* flatten(PTR p) {
    typedef typename deduce<PTR>::type T;
    return reinterpret_cast<T*>(p);
  }

  //====================================================================//

  template <class PTR>
  static void allocate(size_t nelem, PTR* pp) {
    typedef typename deduce<PTR>::type T;
    constexpr size_t N = deduce<PTR>::N;
    Allocator<T> a;
    a.allocate_bytes(reinterpret_cast<void**>(pp), sizeof(T) * nelem * N);
  }

  template <class PTR>
  static void deallocate(PTR p) {
    typedef typename deduce<PTR>::type T;
    constexpr size_t N = deduce<PTR>::N;
    Allocator<T> a;
    a.deallocate_bytes(flatten(p));
  }

  template <class PTR, class... PTRS>
  static void allocate(size_t nelem, PTR* pp, PTRS... pps) {
    allocate(nelem, pp);
    allocate(nelem, pps...);
  }

  template <class PTR, class... PTRS>
  static void deallocate(PTR p, PTRS... ps) {
    deallocate(p);
    deallocate(ps...);
  }

  //====================================================================//

  template <class PTR>
  static void zero(size_t nelem, PTR p) {
    typedef typename deduce<PTR>::type T;
    constexpr size_t N = deduce<PTR>::N;
    Allocator<T> a;
    a.zero_bytes(flatten(p), sizeof(T) * nelem * N);
  }

  template <class PTR, class... PTRS>
  static void zero(size_t nelem, PTR p, PTRS... ps) {
    zero(nelem, p);
    zero(nelem, ps...);
  }

  //====================================================================//

  template <class PTR, class U>
  static void copyin(size_t nelem, PTR dst, const U* src) {
    typedef typename deduce<PTR>::type T;
    constexpr size_t N = deduce<PTR>::N;
    Allocator<T> a;
    a.copyin_array(flatten(dst), flatten(src), nelem * N);
  }

  template <class U, class PTR>
  static void copyout(size_t nelem, U* dst, const PTR src) {
    typedef typename deduce<PTR>::type T;
    constexpr size_t N = deduce<PTR>::N;
    Allocator<T> a;
    a.copyout_array(flatten(dst), flatten(src), nelem * N);
  }

  template <class PTR, class U>
  static void copy(size_t nelem, PTR dst, const U* src) {
    typedef typename deduce<PTR>::type T;
    constexpr size_t N = deduce<PTR>::N;
    Allocator<T> a;
    a.copy_array(flatten(dst), flatten(src), nelem * N);
  }

  //====================================================================//

  template <class PTR, class PTR2>
  static typename deduce<PTR>::type dot(size_t nelem, const PTR ptr,
                                        const PTR2 b) {
    typedef typename deduce<PTR>::type T;
    constexpr size_t N = deduce<PTR>::N;
    typedef typename deduce<PTR2>::type T2;
    static_assert(std::is_same<T, T2>::value, "");
    return mathfunc_detail::dotprod(flatten(ptr), flatten(b), nelem * N);
  }

  //====================================================================//

  template <class FLT, class PTR>
  static void scale(size_t nelem, FLT scal, PTR ptr) {
    typedef typename deduce<PTR>::type T;
    constexpr size_t N = deduce<PTR>::N;
    mathfunc_detail::scale_array(flatten(ptr), scal, nelem * N);
  }

  template <class FLT, class PTR, class... PTRS>
  static void scale(size_t nelem, FLT scal, PTR ptr, PTRS... ptrs) {
    scale(nelem, scal, ptr);
    scale(nelem, scal, ptrs...);
  }
};

template <class T>
using DeviceArrayAllocator = DeviceAllocator<T, DeviceMemory>;
typedef DeviceArray<DeviceArrayAllocator> device_array;

template <class T, size_t N = 1>
using device_pointer = typename device_array::ptr<T, N>::type;
TINKER_NAMESPACE_END

#endif
