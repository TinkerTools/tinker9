#include "compilers.h"
#include "io_print.h"


TINKER_NAMESPACE_BEGIN
std::string cxx_compiler_name()
{
   std::string n = "unknown";
#if defined(__INTEL_COMPILER)
   n = format("icpc %d.%d", __INTEL_COMPILER, __INTEL_COMPILER_BUILD_DATE);

#elif defined(__clang__)
   // xcode clang is different from llvm clang
#   ifdef __apple_build_version__
   n = format("clang++ %d.%d.%d (xcode)", __clang_major__, __clang_minor__,
              __clang_patchlevel__);
#   else
   n = format("clang++ %d.%d.%d (llvm)", __clang_major__, __clang_minor__,
              __clang_patchlevel__);
#   endif

#elif defined(__GNUC__)
   n = format("g++ %d.%d.%d", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);

#endif
   return n;
}
TINKER_NAMESPACE_END
