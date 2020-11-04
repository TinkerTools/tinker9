#pragma once
#if TINKER_CUDART
#   include "macro.h"
#   include <cublas_v2.h>
#   include <cuda_runtime.h>


namespace tinker {
namespace g {
TINKER_EXTERN cudaStream_t s0;
TINKER_EXTERN cudaStream_t s1;
TINKER_EXTERN cudaStream_t spme;
TINKER_EXTERN cublasHandle_t h0;
TINKER_EXTERN cublasHandle_t h1;
}
TINKER_EXTERN void* pinned_buf;
TINKER_EXTERN void* dptr_buf;
TINKER_EXTERN cudaEvent_t pme_event_finish;
}
#endif
