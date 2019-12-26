#pragma once
#include "rc_man.h"
#if TINKER_CUDART
#   include <cublas_v2.h>
#   include <cuda_runtime.h>
#endif


TINKER_NAMESPACE_BEGIN
#if TINKER_CUDART
extern cudaStream_t nonblk;
extern cublasHandle_t h_cublas;
extern cublasHandle_t h_cublas_nonblk;
extern real* pinned_real64;
extern real* dptr_real64;
#endif


void cudalib_data(rc_op);
TINKER_NAMESPACE_END
