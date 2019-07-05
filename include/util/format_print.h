#ifndef TINKER_UTIL_FORMAT_PRINT_H_
#define TINKER_UTIL_FORMAT_PRINT_H_

#include "cxx.h"

// fmtlib
#include <ext/fmt/ostream.h>

TINKER_NAMESPACE_BEGIN
template <class Out, class Fmt, class... Ts>
void print(Out& out, const Fmt& fmtstr, const Ts&... args) {
  fmt::print(out, fmtstr, args...);
}

template <class Fmt, class... Ts>
std::string format(const Fmt& fmtstr, const Ts&... args) {
  return fmt::format(fmtstr, args...);
}
TINKER_NAMESPACE_END

#endif
