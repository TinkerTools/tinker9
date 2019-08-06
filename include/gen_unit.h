#ifndef TINKER_GEN_UNIT_H_
#define TINKER_GEN_UNIT_H_

#include "macro.h"
#include <cassert>
#include <memory>
#include <vector>

TINKER_NAMESPACE_BEGIN
enum class GenericUnitVersion { V0, V1 };

template <GenericUnitVersion VERS>
struct GenericUnitAlloc;

template <>
struct GenericUnitAlloc<GenericUnitVersion::V0> {
  struct Dealloc {
    void operator()(void*) {}
  };

  struct Alloc {
    void operator()(void**, size_t) {}
  };

  struct Copyin {
    void operator()(void*, const void*, size_t) {}
  };
};

/**
 * @brief
 * resource handle
 *
 * analogous to to Fortran I/O unit that can be used as signed integers
 *
 * @tparam VERSION
 * mainly used for identifying whether to allocate memory on device and store
 * the device pointer that corresponds to the host object; can be extended to
 * specify different de/allocation methods
 */
template <class T, GenericUnitVersion VERSION = GenericUnitVersion::V0>
class GenericUnit {
private:
  int m_unit;

  static constexpr int USE_DPTR = (VERSION == GenericUnitVersion::V0 ? 0 : 1);

  typedef std::vector<std::unique_ptr<T>> hostptr_vec;
  static hostptr_vec& hostptrs() {
    static hostptr_vec o;
    return o;
  }

  typedef typename GenericUnitAlloc<VERSION>::Dealloc Dealloc;
  typedef std::vector<std::unique_ptr<T, Dealloc>> dptr_vec;
  static dptr_vec& deviceptrs() {
    assert(USE_DPTR);
    static dptr_vec o;
    return o;
  }

  typedef typename GenericUnitAlloc<VERSION>::Alloc Alloc;
  typedef typename GenericUnitAlloc<VERSION>::Copyin Copyin;

public:
  static int size() {
    if_constexpr(USE_DPTR) assert(hostptrs().size() == deviceptrs().size());
    return hostptrs().size();
  }

  static void clear() {
    // call ~T() here
    hostptrs().clear();
    // call Dealloc(T*) here
    if_constexpr(USE_DPTR) deviceptrs().clear();
  }

  static void resize(int s) {
    assert(!USE_DPTR);
    for (int i = size(); i < s; ++i) {
      hostptrs().emplace_back(new T);
    }
  }

  static GenericUnit alloc_new() {
    hostptrs().emplace_back(new T);
    if_constexpr(USE_DPTR) {
      T* ptr;
      Alloc alloc;
      alloc(reinterpret_cast<void**>(&ptr), sizeof(T));
      Dealloc dealloc;
      deviceptrs().emplace_back(ptr, dealloc);
    }
    return size() - 1;
  }

public:
  GenericUnit()
      : m_unit(-1) {}

  GenericUnit(int u)
      : m_unit(u) {}

  operator int() const { return m_unit; }

  const T& obj() const {
    assert(0 <= m_unit && m_unit < hostptrs().size() &&
           "const T& GenericUnit::obj() const");
    return *hostptrs()[m_unit];
  }

  T& obj() {
    assert(0 <= m_unit && m_unit < hostptrs().size() &&
           "T& GenericUnit::obj()");
    return *hostptrs()[m_unit];
  }

  const T* deviceptr() const {
    assert(0 <= m_unit && m_unit < deviceptrs().size() &&
           "const T* GenericUnit::deviceptr() const");
    return deviceptrs()[m_unit].get();
  }

  T* deviceptr() {
    assert(0 <= m_unit && m_unit < deviceptrs().size() &&
           "T* GenericUnit::deviceptr()");
    return deviceptrs()[m_unit].get();
  }

  void init_deviceptr(const T& hobj) {
    assert(&hobj == &this->obj());
    Copyin copyin;
    copyin(this->deviceptr(), &this->obj(), sizeof(T));
  }
};
TINKER_NAMESPACE_END

#endif
