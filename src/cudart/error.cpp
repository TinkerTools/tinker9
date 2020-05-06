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


namespace tinker {
template <>
std::string translate_error_code<cudaError_t>(cudaError_t error_num)
{
   return std::string(cudaGetErrorString(error_num));
}


template <>
std::string translate_error_code<cublasStatus_t>(cublasStatus_t error_num)
{
   return "";
}


template <>
std::string translate_error_code<cufftResult_t>(cufftResult_t error_num)
{
   return "";
}


template <>
std::string translate_error_code<nvmlReturn_enum>(nvmlReturn_enum error_num)
{
   return "";
}
}
