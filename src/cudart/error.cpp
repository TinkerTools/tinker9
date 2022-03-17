#include "tool/error.h"
#include <cuda_runtime.h>

extern "C"
{
   struct cublasStatus_t
   {};
   struct cufftResult_t
   {};
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
