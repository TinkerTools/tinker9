#include "tool/compilers.h"
#include "tool/io.h"

namespace tinker {
std::string cudaCompilerName()
{
   return format("nvcc %d.%d.%d", __CUDACC_VER_MAJOR__, __CUDACC_VER_MINOR__, __CUDACC_VER_BUILD__);
}
}
