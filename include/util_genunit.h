#ifndef TINKER_UTIL_GENUNIT_H_
#define TINKER_UTIL_GENUNIT_H_

#include "util_macro.h"
#include <cassert>
#include <cstring>
#include <vector>

TINKER_NAMESPACE_BEGIN
/**
 * @brief
 * can be used as signed integers, similar to Fortran I/O unit
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
    hostobjs_().resize(s, T());
    if_constexpr(USE_DPTR) deviceptrs_().resize(s, nullptr);
  }

  static void emplace_back(T&& v) {
    hostobjs_().emplace_back(std::forward<T>(v));
    if_constexpr(USE_DPTR) deviceptrs_().emplace_back(nullptr);
  }

public:
  GenericUnit()
      : unit_(-1) {}
  GenericUnit(int u)
      : unit_(u) {}
  operator int() const { return unit_; }

  const T& obj() const {
#if TINKER_DEBUG
    return hostobjs_().at(unit_);
#else
    return hostobjs_()[unit_];
#endif
  }

  T& obj() {
#if TINKER_DEBUG
    return hostobjs_().at(unit_);
#else
    return hostobjs_()[unit_];
#endif
  }

  const T* deviceptr() const {
#if TINKER_DEBUG
    return deviceptrs_().at(unit_);
#else
    return deviceptrs_()[unit_];
#endif
  }

  T*& deviceptr() {
#if TINKER_DEBUG
    return deviceptrs_().at(unit_);
#else
    return deviceptrs_()[unit_];
#endif
  }
};
TINKER_NAMESPACE_END

#endif
