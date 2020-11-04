#pragma once
#if TINKER_CUDART
#   include "macro.h"
#   include <cublas_v2.h>
#   include <cuda_runtime.h>


namespace tinker {
namespace g {
/// \ingroup async
/// \brief CUDA stream for the default OpenACC async queue.
TINKER_EXTERN cudaStream_t s0;
/// \ingroup async
/// \brief CUDA stream for the default OpenACC sync queue.
TINKER_EXTERN cudaStream_t s1;
/// \ingroup async
/// \brief CUDA stream for the OpenACC async %PME queue.
TINKER_EXTERN cudaStream_t spme;
/// \ingroup async
/// \brief CUDA BLAS handle using #s0.
TINKER_EXTERN cublasHandle_t h0;
/// \ingroup async
/// \brief CUDA BLAS handle using #s1.
TINKER_EXTERN cublasHandle_t h1;
}
TINKER_EXTERN void* pinned_buf;
TINKER_EXTERN void* dptr_buf;
TINKER_EXTERN cudaEvent_t pme_event_finish;
}
#endif
