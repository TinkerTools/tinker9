#if defined(__APPLE__) || defined(__linux__)
#  include "util_io.h"
#  include <cxxabi.h>
#  include <execinfo.h>

TINKER_NAMESPACE_BEGIN
namespace detail_ {
enum class BackTraceOS { macOS, Linux };

template <BackTraceOS os>
void print_backtrace_tmpl(std::ostream& fp) {
  const int max_frames = 128;
  void* callstack[max_frames];
  int frames = backtrace(callstack, max_frames);
  char** strs = backtrace_symbols(callstack, frames);

  char* demangled_name;
  int status;
  std::string num, caller, callee;
  const char* f1 = " Backtrace\n";
  const char* f2 = " {:>4s}  {:60.60s}  {:s}\n";
  const int tolerance = 20;
  print(fp, f1);

  for (auto i = 1; i < frames; ++i) {
    std::string tmp = strs[i];

    // e.g.
    // 3   libdyld.dylib                       0x00007fffbc358235 start + 1
    // [0] [1]                                 [2]                [3]     [5]
    if_constexpr(os == BackTraceOS::macOS) {
      auto vs = Text::split(tmp);
      num = vs.at(0);
      caller = vs.at(1);
      callee = vs.at(3);
    }

    // e.g.
    // /lib/x86_64-linux-gnu/libc.so.6(__libc_start_main+0xe7) [0x7f4c7dd9db97]
    // [0]                             [1]               [2]    [3]
    else if_constexpr(os == BackTraceOS::Linux) {
      Text::replace(tmp, "()[]+", ' ');
      auto vs = Text::split(tmp);
      num = std::to_string(i);
      caller = vs.at(0);
      callee = vs.at(1);
    }

    demangled_name = abi::__cxa_demangle(callee.c_str(), 0, 0, &status);
    if (!status) {
      // This name CAN be demangled.
      std::string tmp =
          ((tolerance + callee.length()) >= std::strlen(demangled_name))
          ? demangled_name
          : callee;
      callee = tmp;
      free(demangled_name);
    }
    print(fp, f2, num, caller, callee);
  }
  free(strs);

  fp << std::flush;
}
}

void print_backtrace(std::ostream& fp) {
#  if defined(__APPLE__)
  detail_::print_backtrace_tmpl<detail_::BackTraceOS::macOS>(fp);
#  elif defined(__linux__)
  detail_::print_backtrace_tmpl<detail_::BackTraceOS::Linux>(fp);
#  endif
}
TINKER_NAMESPACE_END

#endif
