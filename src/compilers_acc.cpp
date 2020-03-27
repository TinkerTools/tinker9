#include "compilers.h"
#include "io_print.h"


TINKER_NAMESPACE_BEGIN
#if TINKER_CUDART
std::string acc_compiler_name()
{
   return format("pgc++ %d.%d.%d", __PGIC__, __PGIC_MINOR__,
                 __PGIC_PATCHLEVEL__);
}
#endif
TINKER_NAMESPACE_END
