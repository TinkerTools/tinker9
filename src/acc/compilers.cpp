#include "tool/compilers.h"
#include "tool/ioprint.h"

namespace tinker {
std::string accCompilerName()
{
#if TINKER_GPULANG_OPENACC
   if (__PGIC__ <= 19)
      return format("pgc++ %d.%d.%d", __PGIC__, __PGIC_MINOR__, __PGIC_PATCHLEVEL__);
   else
      return format("nvc++ %d.%d.%d", __PGIC__, __PGIC_MINOR__, __PGIC_PATCHLEVEL__);
#else
   return "Unused";
#endif
}
}
