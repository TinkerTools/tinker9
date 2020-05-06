#include "compilers.h"
#include "io_print.h"


namespace tinker {
std::string cuda_compiler_name()
{
   return format("nvcc %d.%d.%d", __CUDACC_VER_MAJOR__, __CUDACC_VER_MINOR__,
                 __CUDACC_VER_BUILD__);
}
}
