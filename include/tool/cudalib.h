#pragma once
#include "tool/macro.h"
#include "tool/rcman.h"

namespace tinker {
/// \ingroup platform
/// \brief Sets up the CUDA variables including but not limited to CUDA streams,
/// CUDA library handles, CUDA memory buffers, and integer units for the OpenACC queues.
void cudalibData(RcOp);
}

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

#if TINKER_CUDART
#include <cublas_v2.h>
#include <cuda_runtime.h>

namespace tinker {
namespace g {
TINKER_EXTERN cudaStream_t s0;   ///< CUDA stream for the default OpenACC async queue. \ingroup async
TINKER_EXTERN cudaStream_t s1;   ///< CUDA stream for the default OpenACC sync queue. \ingroup async
TINKER_EXTERN cudaStream_t spme; ///< CUDA stream for the OpenACC async %PME queue. \ingroup pme
TINKER_EXTERN cublasHandle_t h0; ///< CUDA BLAS handle using #s0. \ingroup async
TINKER_EXTERN cublasHandle_t h1; ///< CUDA BLAS handle using #s1. \ingroup async
}

TINKER_EXTERN void* pinned_buf; ///< Preallocated pinned host memory. \ingroup platform
TINKER_EXTERN void* dptr_buf;   ///< Preallocated device memory. \ingroup platform

TINKER_EXTERN cudaEvent_t pme_event_start;  ///< \ingroup pme
TINKER_EXTERN cudaEvent_t pme_event_finish; ///< \ingroup pme
}
#endif
