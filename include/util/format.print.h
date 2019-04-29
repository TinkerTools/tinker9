#ifndef TINKER_UTIL_FORMAT_PRINT_H_
#define TINKER_UTIL_FORMAT_PRINT_H_

#include "cxx.h"

// fmtlib
#include <ext/fmt/ostream.h>

TINKER_NAMESPACE_BEGIN
template <class __Out, class __Fmt, class... __Ts>
void print(__Out& _out, const __Fmt& _fmtstr, const __Ts&... _args) {
  fmt::print(_out, _fmtstr, _args...);
}

template <class __Fmt, class... __Ts>
std::string format(const __Fmt& _fmtstr, const __Ts&... _args) {
  return fmt::format(_fmtstr, _args...);
}
TINKER_NAMESPACE_END

#endif
