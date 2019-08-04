#ifndef TINKER_UTIL_GENUNIT_H_
#define TINKER_UTIL_GENUNIT_H_

#include "util_macro.h"
#include <vector>

TINKER_NAMESPACE_BEGIN
template <class T>
class GenericUnit {
private:
  int unit_;

public:
  GenericUnit()
      : unit_(-1) {}
  GenericUnit(int u)
      : unit_(u) {}
  bool operator<(const GenericUnit<T>& u) const { return unit_ < u.unit_; }
  bool operator==(const GenericUnit<T>& u) const { return unit_ == u.unit_; }
  int unit() const { return unit_; }

public:
  static std::vector<T>& all_objs() {
    static std::vector<T> o;
    return o;
  }

  static std::vector<T*>& all_deviceptrs() {
    static std::vector<T*> o;
    return o;
  }
};
TINKER_NAMESPACE_END

#endif
