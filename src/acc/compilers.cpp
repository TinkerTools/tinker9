#include "tool/compilers.h"
#include "tool/io.h"

namespace tinker {
#if TINKER_CUDART
std::string accCompilerName()
{
   if (__PGIC__ <= 19)
      return format("pgc++ %d.%d.%d", __PGIC__, __PGIC_MINOR__, __PGIC_PATCHLEVEL__);
   else
      return format("nvc++ %d.%d.%d", __PGIC__, __PGIC_MINOR__, __PGIC_PATCHLEVEL__);
}
#endif
}
