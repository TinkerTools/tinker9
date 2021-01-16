#include "tool/compilers.h"
#include "tool/io_print.h"


namespace tinker {
#if TINKER_CUDART
std::string acc_compiler_name()
{
   if (__PGIC__ <= 19)
      return format("pgc++ %d.%d.%d", __PGIC__, __PGIC_MINOR__,
                    __PGIC_PATCHLEVEL__);
   else
      return format("nvc++ %d.%d.%d", __PGIC__, __PGIC_MINOR__,
                    __PGIC_PATCHLEVEL__);
}
#endif
}
