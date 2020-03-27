#include "compilers.h"
#include "io_print.h"


TINKER_NAMESPACE_BEGIN
std::string cuda_compiler_name()
{
   return format("nvcc %d.%d.%d", __CUDACC_VER_MAJOR__, __CUDACC_VER_MINOR__,
                 __CUDACC_VER_BUILD__);
}
TINKER_NAMESPACE_END
