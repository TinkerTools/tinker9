#pragma once
#if TINKER_CUDART
#   include "macro.h"
#   include <cublas_v2.h>
#   include <cuda_runtime.h>


namespace tinker {
TINKER_EXTERN cudaStream_t nonblk;
TINKER_EXTERN cublasHandle_t h_cublas;
TINKER_EXTERN cublasHandle_t h_cublas_nonblk;
TINKER_EXTERN void* pinned_buf;
TINKER_EXTERN void* dptr_buf;
}
#endif
