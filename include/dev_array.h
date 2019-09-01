#ifndef TINKER_DEV_ARRAY_H_
#define TINKER_DEV_ARRAY_H_

#include "dev_mem.h"

TINKER_NAMESPACE_BEGIN
class device_array {
private:
  template <class T>
  using allocator = DeviceAllocator<T>;

  template <class PTR>
  struct deduce;

  template <class T>
  struct deduce<T*> {
    typedef T type;
    static constexpr size_t N = 1;
  };

  template <class T, size_t M>
  struct deduce<T (*)[M]> {
    typedef T type;
    static constexpr size_t N = M;
  };

public:
  template <class T, size_t N = 1>
  struct ptr;

  template <class T>
  struct ptr<T> {
    typedef T* type;
    static T* flatten(type p) { return &p[0]; }
  };

  template <class T, size_t N>
  struct ptr {
    static_assert(N > 1, "");
    typedef T (*type)[N];
    static T* flatten(type p) { return &p[0][0]; }
  };

  //====================================================================//

  template <class PTR>
  static void allocate(PTR* pp, size_t nelem) {
    typedef typename deduce<PTR>::type T;
    constexpr size_t N = deduce<PTR>::N;
    allocator<T> a;
    a.allocate_bytes(reinterpret_cast<void**>(pp), sizeof(T) * nelem * N);
  }

  template <class PTR>
  static void deallocate(PTR p) {
    typedef typename deduce<PTR>::type T;
    constexpr size_t N = deduce<PTR>::N;
    allocator<T> a;
    a.deallocate_bytes(ptr<T, N>::flatten(p));
  }

  template <class PTR, class... PTRS>
  static void deallocate(PTR p, PTRS... ps) {
    deallocate(p);
    deallocate(ps...);
  }

  //====================================================================//

  template <class PTR>
  static void zero(PTR p, size_t nelem) {
    typedef typename deduce<PTR>::type T;
    constexpr size_t N = deduce<PTR>::N;
    allocator<T> a;
    a.zero_bytes(ptr<T, N>::flatten(p), sizeof(T) * nelem * N);
  }

  //====================================================================//

  template <class PTR, class U>
  static void copyin(PTR dst, const U* src, size_t nelem) {
    typedef typename deduce<PTR>::type T;
    constexpr size_t N = deduce<PTR>::N;
    allocator<T> a;
    a.copyin_array(ptr<T, N>::flatten(dst), src, nelem * N);
  }
};
TINKER_NAMESPACE_END

#endif
