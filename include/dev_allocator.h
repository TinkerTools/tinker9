#ifndef TINKER_DEV_ALLOCATOR_H_
#define TINKER_DEV_ALLOCATOR_H_

#include "dev_memory.h"
#include <stdexcept>

TINKER_NAMESPACE_BEGIN
template <class T, class ByteOp>
struct DeviceAllocator : public ByteOp {
  typedef T value_type;

  ~DeviceAllocator() = default;

  DeviceAllocator() = default;

  template <class U>
  constexpr DeviceAllocator(const DeviceAllocator<U, ByteOp>&) noexcept {}

  void deallocate(T* p, size_t) noexcept { ByteOp::deallocate_bytes(p); }

  T* allocate(size_t nelem) {
    void* p = nullptr;

    if (nelem > size_t(-1) / sizeof(T))
      goto bad_flag;

    ByteOp::allocate_bytes(&p, nelem * sizeof(T));
    if (p)
      return reinterpret_cast<T*>(p);
    else
    bad_flag:
      throw std::bad_alloc();
  }
};

template <class T1, class T2, class ByteOp>
bool operator==(const DeviceAllocator<T1, ByteOp>&,
                const DeviceAllocator<T2, ByteOp>&) {
  return true;
}

template <class T1, class T2, class ByteOp>
bool operator!=(const DeviceAllocator<T1, ByteOp>&,
                const DeviceAllocator<T2, ByteOp>&) {
  return false;
}
TINKER_NAMESPACE_END

#endif
