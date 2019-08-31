#ifndef TINKER_DEVICE_MEMORY_H_
#define TINKER_DEVICE_MEMORY_H_

#include "macro.h"
#include <cstring>
#include <stdexcept>
#include <type_traits>
#include <vector>

TINKER_NAMESPACE_BEGIN
struct DeviceMemory {
  static void copyin_bytes(void* dst, const void* src, size_t nbytes);

  struct CopyIn {
    void operator()(void* dst, const void* src, size_t nbytes) {
      copyin_bytes(dst, src, nbytes);
    }
  };

  template <class DT, class ST>
  static void copyin_array(DT* dst, const ST* src, size_t nelem) {
    constexpr size_t ds = sizeof(DT); // device type
    constexpr size_t ss = sizeof(ST); // host type
    static_assert(ds <= ss, "invalid if sizeof(DstType) > sizeof(SrcType)");

    size_t size = ds * nelem;
    if_constexpr(ds == ss) { copyin_bytes(dst, src, size); }
    else if_constexpr(ds < ss) {
      std::vector<DT> buf(nelem);
      for (size_t i = 0; i < nelem; ++i)
        buf[i] = src[i];
      copyin_bytes(dst, buf.data(), size);
    }
  }

  static void copyout_bytes(void* dst, const void* src, size_t nbytes);

  struct CopyOut {
    void operator()(void* dst, const void* src, size_t nbytes) {
      copyout_bytes(dst, src, nbytes);
    }
  };

  template <class DT, class ST>
  static void copyout_array(DT* dst, const ST* src, size_t nelem) {
    constexpr size_t ds = sizeof(DT); // host type
    constexpr size_t ss = sizeof(ST); // device type
    static_assert(ds >= ss, "invalid if sizeof(SrcType) > sizeof(DstType)");

    size_t size = ss * nelem;
    if_constexpr(ds == ss) { copyout_bytes(dst, src, size); }
    else if_constexpr(ds > ss) {
      std::vector<ST> buf(nelem);
      copyout_bytes(buf.data(), src, size);
      for (size_t i = 0; i < nelem; ++i)
        dst[i] = buf[i];
    }
  }

  static void copy_bytes(void* dst, const void* src, size_t nbytes);

  struct Copy {
    void operator()(void* dst, const void* src, size_t nbytes) {
      copy_bytes(dst, src, nbytes);
    }
  };

  template <class DT, class ST>
  static void copy_array(DT* dst, const ST* src, size_t nelem) {
    static_assert(std::is_same<DT, ST>::value, "");
    size_t size = sizeof(ST) * nelem;
    copy_bytes(dst, src, size);
  }

  static void zero_bytes(void* ptr, size_t nbytes);

  struct Zero {
    void operator()(void* ptr, size_t nbytes) { zero_bytes(ptr, nbytes); }
  };

  template <class T>
  static void zero_array(T* dst, size_t nelem) {
    size_t size = sizeof(T) * nelem;
    zero_bytes(dst, size);
  }

  static void deallocate_bytes(void* ptr);

  struct Deallocate {
    void operator()(void* ptr) { deallocate_bytes(ptr); }
  };

  static void allocate_bytes(void** pptr, size_t nbytes);

  template <
      class T,
      class = typename std::enable_if<!std::is_same<T, void>::value>::type>
  void allocate_types(T** pptr, size_t nbytes) {
    return allocate_bytes(reinterpret_cast<void**>(pptr), nbytes);
  }

  struct Allocate {
    void operator()(void** pptr, size_t nbytes) {
      allocate_bytes(pptr, nbytes);
    }
  };
};

template <class T>
struct DeviceAllocator : public DeviceMemory {
  typedef T value_type;

  ~DeviceAllocator() = default;

  DeviceAllocator() = default;

  template <class U>
  constexpr DeviceAllocator(const DeviceAllocator<U>&) noexcept {}

  void deallocate(T* p, size_t) noexcept { deallocate_bytes(p); }

  T* allocate(size_t nelem) {
    void* p = nullptr;

    if (nelem > size_t(-1) / sizeof(T))
      goto bad_flag;

    allocate_bytes(&p, nelem * sizeof(T));
    if (p)
      return reinterpret_cast<T*>(p);
    else
    bad_flag:
      throw std::bad_alloc();
  }
};

template <class T1, class T2>
bool operator==(const DeviceAllocator<T1>&, const DeviceAllocator<T2>&) {
  return true;
}

template <class T1, class T2>
bool operator!=(const DeviceAllocator<T1>&, const DeviceAllocator<T2>&) {
  return false;
}
TINKER_NAMESPACE_END

#endif
