#ifndef TINKER_UTIL_MATH_H_
#define TINKER_UTIL_MATH_H_

#include "macro.h"
#include <algorithm>

TINKER_NAMESPACE_BEGIN
template <class T>
T max_of(T a) {
  return a;
}

template <class T>
T max_of(T a, T b) {
  return (a < b) ? b : a;
}

template <class T, class... Ts>
T max_of(T a, T b, Ts... cs) {
  return max_of(max_of(a, b), cs...);
}

template <class T>
T min_of(T a) {
  return a;
}

template <class T>
T min_of(T a, T b) {
  return (a < b) ? a : b;
}

template <class T, class... Ts>
T min_of(T a, T b, Ts... cs) {
  return min_of(min_of(a, b), cs...);
}
TINKER_NAMESPACE_END

#endif
