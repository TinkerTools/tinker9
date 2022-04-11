#include "tool/compilers.h"
#include "tool/ioprint.h"

namespace tinker {
#if TINKER_GPULANG_OPENACC
std::string accCompilerName()
{
   if (__PGIC__ <= 19)
      return format("pgc++ %d.%d.%d", __PGIC__, __PGIC_MINOR__, __PGIC_PATCHLEVEL__);
   else
      return format("nvc++ %d.%d.%d", __PGIC__, __PGIC_MINOR__, __PGIC_PATCHLEVEL__);
}
#endif
}
