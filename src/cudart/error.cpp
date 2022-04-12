#include "tool/error.h"
#include "tool/macro.h"
#include <cuda_runtime.h>

TINKER_DECL_EXTN("C")
{
   struct cublasStatus_t
   {
      int foo;
   };

   struct cufftResult_t
   {
      int foo;
   };
}

namespace tinker {
template <>
std::string translateErrorCode<cudaError_t>(cudaError_t error_num)
{
   return std::string(cudaGetErrorString(error_num));
}

template <>
std::string translateErrorCode<cublasStatus_t>(cublasStatus_t error_num)
{
   return "";
}

template <>
std::string translateErrorCode<cufftResult_t>(cufftResult_t error_num)
{
   return "";
}
}
