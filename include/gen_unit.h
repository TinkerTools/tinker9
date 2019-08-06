#ifndef TINKER_GEN_UNIT_H_
#define TINKER_GEN_UNIT_H_

#include "macro.h"
#include <cassert>
#include <cstring>
#include <vector>

TINKER_NAMESPACE_BEGIN
/**
 * @brief
 * resource handle: can be used as signed integers, similar to Fortran I/O unit
 *
 * @tparam USE_DPTR
 * whether to allocate memory on device and store the device pointer that
 * corresponds to the host object
 */
template <class T, int USE_DPTR = 0>
class GenericUnit {
private:
  int unit_;

  static std::vector<T>& hostobjs_() {
    static std::vector<T> o;
    return o;
  }

  static std::vector<T*>& deviceptrs_() {
    assert(USE_DPTR);
    static std::vector<T*> o;
    return o;
  }

public:
  static size_t size() {
    if_constexpr(USE_DPTR) assert(hostobjs_().size() == deviceptrs_().size());
    return hostobjs_().size();
  }

  static void clear() {
    hostobjs_().clear();
    if_constexpr(USE_DPTR) deviceptrs_().clear();
  }

  static void resize(size_t s) {
    hostobjs_().resize(s);
    if_constexpr(USE_DPTR) deviceptrs_().resize(s);
  }

  static GenericUnit add_new() {
    hostobjs_().emplace_back(T());
    if_constexpr(USE_DPTR) deviceptrs_().emplace_back(nullptr);
    return size() - 1;
  }

public:
  GenericUnit()
      : unit_(-1) {}
  GenericUnit(int u)
      : unit_(u) {}
  operator int() const { return unit_; }

  const T& obj() const {
#if TINKER_DEBUG
    assert(0 <= unit_ && unit_ < hostobjs_().size() &&
           "GenericUnit::obj() const");
#endif
    return hostobjs_()[unit_];
  }

  T& obj() {
#if TINKER_DEBUG
    assert(0 <= unit_ && unit_ < hostobjs_().size() && "GenericUnit::obj()");
#endif
    return hostobjs_()[unit_];
  }

  const T* deviceptr() const {
#if TINKER_DEBUG
    assert(0 <= unit_ && unit_ < deviceptrs_().size() &&
           "GenericUnit::deviceptr() const");
#endif
    return deviceptrs_()[unit_];
  }

  T*& deviceptr() {
#if TINKER_DEBUG
    assert(0 <= unit_ && unit_ < deviceptrs_().size() &&
           "GenericUnit::deviceptr()");
#endif
    return deviceptrs_()[unit_];
  }
};
TINKER_NAMESPACE_END

#endif
