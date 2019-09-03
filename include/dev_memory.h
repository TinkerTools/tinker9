#ifndef TINKER_DEV_MEMORY_H_
#define TINKER_DEV_MEMORY_H_

#include "macro.h"
#include <cstring>
#include <type_traits>
#include <vector>

TINKER_NAMESPACE_BEGIN
void device_memory_copyin_bytes(void* dst, const void* src, size_t nbytes);
void device_memory_copyout_bytes(void* dst, const void* src, size_t nbytes);
void device_memory_copy_bytes(void* dst, const void* src, size_t nbytes);
void device_memory_zero_bytes(void* dst, size_t nbytes);
void device_memory_deallocate_bytes(void* ptr);
void device_memory_allocate_bytes(void** pptr, size_t nbytes);

struct DeviceMemory {
  static void copyin_bytes(void* dst, const void* src, size_t nbytes) {
    device_memory_copyin_bytes(dst, src, nbytes);
  }

  static void copyout_bytes(void* dst, const void* src, size_t nbytes) {
    device_memory_copyout_bytes(dst, src, nbytes);
  }

  static void copy_bytes(void* dst, const void* src, size_t nbytes) {
    device_memory_copyout_bytes(dst, src, nbytes);
  }

  static void zero_bytes(void* dst, size_t nbytes) {
    device_memory_zero_bytes(dst, nbytes);
  }

  static void deallocate_bytes(void* ptr) {
    device_memory_deallocate_bytes(ptr);
  }

  static void allocate_bytes(void** pptr, size_t nbytes) {
    device_memory_allocate_bytes(pptr, nbytes);
  }

  //====================================================================//

  template <
      class T,
      class = typename std::enable_if<!std::is_same<T, void>::value>::type>
  static void allocate_bytes(T** pptr, size_t nbytes) {
    return allocate_bytes(reinterpret_cast<void**>(pptr), nbytes);
  }

  //====================================================================//

  struct CopyIn {
    void operator()(void* dst, const void* src, size_t nbytes) {
      copyin_bytes(dst, src, nbytes);
    }
  };

  struct CopyOut {
    void operator()(void* dst, const void* src, size_t nbytes) {
      copyout_bytes(dst, src, nbytes);
    }
  };

  struct Copy {
    void operator()(void* dst, const void* src, size_t nbytes) {
      copy_bytes(dst, src, nbytes);
    }
  };

  struct Zero {
    void operator()(void* ptr, size_t nbytes) { zero_bytes(ptr, nbytes); }
  };

  struct Deallocate {
    void operator()(void* ptr) { deallocate_bytes(ptr); }
  };

  struct Allocate {
    void operator()(void** pptr, size_t nbytes) {
      allocate_bytes(pptr, nbytes);
    }
  };

  //====================================================================//

  template <class T>
  static void check_type() {
    static_assert(std::is_enum<T>::value || std::is_integral<T>::value ||
                      std::is_floating_point<T>::value ||
                      std::is_trivial<T>::value,
                  "");
  }

  //====================================================================//

  template <class DT, class ST>
  static void copyin_array(DT* dst, const ST* src, size_t nelem) {
    check_type<DT>();
    check_type<ST>();
    constexpr size_t ds = sizeof(DT); // device type
    constexpr size_t ss = sizeof(ST); // host type

    size_t size = ds * nelem;
    if_constexpr(ds == ss) { copyin_bytes(dst, src, size); }
    else {
      std::vector<DT> buf(nelem);
      for (size_t i = 0; i < nelem; ++i)
        buf[i] = src[i];
      copyin_bytes(dst, buf.data(), size);
    }
  }

  template <class DT, class ST>
  static void copyout_array(DT* dst, const ST* src, size_t nelem) {
    check_type<DT>();
    check_type<ST>();
    constexpr size_t ds = sizeof(DT); // host type
    constexpr size_t ss = sizeof(ST); // device type

    size_t size = ss * nelem;
    if_constexpr(ds == ss) { copyout_bytes(dst, src, size); }
    else {
      std::vector<ST> buf(nelem);
      copyout_bytes(buf.data(), src, size);
      for (size_t i = 0; i < nelem; ++i)
        dst[i] = buf[i];
    }
  }

  template <class DT, class ST>
  static void copy_array(DT* dst, const ST* src, size_t nelem) {
    check_type<DT>();
    check_type<ST>();
    static_assert(std::is_same<DT, ST>::value, "");
    size_t size = sizeof(ST) * nelem;
    copy_bytes(dst, src, size);
  }

  template <class T>
  static void zero_array(T* dst, size_t nelem) {
    check_type<T>();
    size_t size = sizeof(T) * nelem;
    zero_bytes(dst, size);
  }

  //====================================================================//

  template <class DT, class ST>
  static void copyin_array2(size_t idx0, size_t ndim, DT* dst, const ST* src,
                            size_t nelem) {
    check_type<DT>();
    check_type<ST>();
    std::vector<DT> buf(nelem);
    for (size_t i = 0; i < nelem; ++i)
      buf[i] = src[ndim * i + idx0];
    copyin_array(dst, buf.data(), nelem);
  }

  template <class DT, class ST>
  static void copyout_array2(size_t idx0, size_t ndim, DT* dst, const ST* src,
                             size_t nelem) {
    check_type<DT>();
    check_type<ST>();
    std::vector<ST> buf(nelem);
    copyout_array(buf.data(), src, nelem);
    for (size_t i = 0; i < nelem; ++i)
      dst[ndim * i + idx0] = buf[i];
  }
};
TINKER_NAMESPACE_END

#endif
