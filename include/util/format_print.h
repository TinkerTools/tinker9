#ifndef TINKER_UTIL_FORMAT_PRINT_H_
#define TINKER_UTIL_FORMAT_PRINT_H_

#include "cxx.h"

// fmtlib
#include <ext/fmt/ostream.h>

TINKER_NAMESPACE_BEGIN
template <class Out_, class Fmt_, class... Ts_>
void print(Out_& out_, const Fmt_& fmtstr_, const Ts_&... args_) {
  fmt::print(out_, fmtstr_, args_...);
}

template <class Fmt_, class... Ts_>
std::string format(const Fmt_& fmtstr_, const Ts_&... args_) {
  return fmt::format(fmtstr_, args_...);
}
TINKER_NAMESPACE_END

#endif
