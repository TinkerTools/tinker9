#include "error.h"
#include <cuda_runtime.h>


extern "C"
{
   struct cublasStatus_t
   {};
   struct cufftResult_t
   {};
   struct nvmlReturn_enum
   {};
}


TINKER_NAMESPACE_BEGIN
template <>
std::string translate_error_message<cudaError_t>(cudaError_t error_num)
{
   return std::string(cudaGetErrorString(error_num));
}


template <>
std::string translate_error_message<cublasStatus_t>(cublasStatus_t error_num)
{
   return "";
}


template <>
std::string translate_error_message<cufftResult_t>(cufftResult_t error_num)
{
   return "";
}


template <>
std::string translate_error_message<nvmlReturn_enum>(nvmlReturn_enum error_num)
{
   return "";
}
TINKER_NAMESPACE_END
