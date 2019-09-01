#ifndef TINKER_GEN_UNIT_H_
#define TINKER_GEN_UNIT_H_

#include "dev_mem.h"
#include <cassert>
#include <memory>
#include <vector>

TINKER_NAMESPACE_BEGIN
enum class GenericUnitVersion { DisableOnDevice, EnableOnDevice };

template <GenericUnitVersion VERS>
struct GenericUnitAlloc;

template <>
struct GenericUnitAlloc<GenericUnitVersion::DisableOnDevice> {
  struct Deallocate {
    void operator()(void*) {}
  };

  struct Allocate {
    void operator()(void**, size_t) {}
  };

  struct CopyIn {
    void operator()(void*, const void*, size_t) {}
  };
};

template <>
struct GenericUnitAlloc<GenericUnitVersion::EnableOnDevice> {
  typedef DeviceMemory::Deallocate Deallocate;
  typedef DeviceMemory::Allocate Allocate;
  typedef DeviceMemory::CopyIn CopyIn;
};

/**
 * @brief
 * resource handle
 *
 * analogous to to fortran i/o unit that can be used as signed integers
 *
 * @tparam VERSION
 * mainly to allow allocate memory on device and store the device pointer that
 * corresponds to the host object; can be extended to specify different
 * de/allocation methods
 */
template <class T,
          GenericUnitVersion VERSION = GenericUnitVersion::DisableOnDevice>
class GenericUnit {
private:
  int unit;

  static constexpr int USE_DPTR =
      (VERSION == GenericUnitVersion::DisableOnDevice ? 0 : 1);

  // The host vector will almost definitely expand its capacity, so if you don't
  // want to implement the move constructors and/or the copy constructors of all
  // the possible type T, don't change vector of host pointers to vector of host
  // objects.
  typedef std::vector<std::unique_ptr<T>> hostptr_vec;
  static hostptr_vec& hostptrs() {
    static hostptr_vec o;
    return o;
  }

  const T& obj() const {
    assert(0 <= unit && unit < hostptrs().size() &&
           "const T& GenericUnit::obj() const");
    return *hostptrs()[unit];
  }

  T& obj() {
    assert(0 <= unit && unit < hostptrs().size() && "T& GenericUnit::obj()");
    return *hostptrs()[unit];
  }

  typedef typename GenericUnitAlloc<VERSION>::Deallocate Deallocate;
  typedef std::vector<std::unique_ptr<T, Deallocate>> deviceptr_vec;
  static deviceptr_vec& deviceptrs() {
    assert(USE_DPTR);
    static deviceptr_vec o;
    return o;
  }

  typedef typename GenericUnitAlloc<VERSION>::Allocate Allocate;
  typedef typename GenericUnitAlloc<VERSION>::CopyIn CopyIn;

public:
  /// @brief
  /// get the number of open units
  static int size() {
    if_constexpr(USE_DPTR) assert(hostptrs().size() == deviceptrs().size());
    return hostptrs().size();
  }

  /// @brief
  /// close all of the units and reset @c size() to 0
  static void clear() {
    // call ~T() here
    hostptrs().clear();
    // call Dealloc(T*) here
    if_constexpr(USE_DPTR) deviceptrs().clear();
  }

  /// @brief
  /// resize the capacity for the objects on host;
  /// cannot be called if device pointers are used
  static void resize(int s) {
    assert(!USE_DPTR);
    for (int i = size(); i < s; ++i)
      hostptrs().emplace_back(new T);
  }

  /// @brief
  /// similar to opening a new fortran i/o unit
  ///
  /// @return
  /// the new unit
  static GenericUnit open() {
    hostptrs().emplace_back(new T);
    if_constexpr(USE_DPTR) {
      T* ptr;
      Allocate alloc;
      alloc(reinterpret_cast<void**>(&ptr), sizeof(T));
      Deallocate dealloc;
      deviceptrs().emplace_back(ptr, dealloc);
    }
    return size() - 1;
  }

public:
  GenericUnit()
      : unit(-1) {}

  GenericUnit(int u)
      : unit(u) {}

  operator int() const { return unit; }

  bool valid() const { return unit >= 0; }
  void close() { unit = -1; }

  /// @brief
  /// get the (const) reference to the object on host
  /// @{
  const T& operator*() const { return obj(); }
  T& operator*() { return obj(); }
  /// @}

  /// @brief
  /// get the (const) pointer to the object on host
  /// @{
  const T* operator->() const { return &obj(); }
  T* operator->() { return &obj(); }
  /// @}

  /// @brief
  /// get device pointer to the object
  /// @{
  const T* deviceptr() const {
    assert(0 <= unit && unit < deviceptrs().size() &&
           "const T* GenericUnit::deviceptr() const");
    return deviceptrs()[unit].get();
  }

  T* deviceptr() {
    assert(0 <= unit && unit < deviceptrs().size() &&
           "T* GenericUnit::deviceptr()");
    return deviceptrs()[unit].get();
  }
  /// @}

  /// @brief
  /// initialized the object on device by an object on host
  ///
  /// @param hobj
  /// the reference to the same object on host
  /// that can be accessed by the same unit number
  void init_deviceptr(const T& hobj) {
    assert(&hobj == &this->obj());
    CopyIn copyin;
    copyin(this->deviceptr(), &this->obj(), sizeof(T));
  }
};
TINKER_NAMESPACE_END

#endif
