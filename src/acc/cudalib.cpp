#include "tool/cudalib.h"
#include "glob.accasync.h"
#include "platform.h"
#if TINKER_CUDART
#   include "tool/error.h"
#   include "tool/gpu_card.h"
#   include <cuda_profiler_api.h>
#   include <openacc.h>
#endif

namespace tinker {
void cudalib_data(rc_op op)
{
#if TINKER_CUDART
   if (op & rc_dealloc) {
      check_rt(cudaProfilerStop());

      g::q0 = -42;
      g::q1 = -42;
      g::s0 = nullptr;
      g::s1 = nullptr;
      check_rt(cublasDestroy(g::h0));
      check_rt(cublasDestroy(g::h1));
      check_rt(cudaFreeHost(pinned_buf));
      check_rt(cudaFree(dptr_buf));

      use_pme_stream = false;
      g::spme = nullptr;
      g::qpme = -42;
      check_rt(cudaEventDestroy(pme_event_start));
      check_rt(cudaEventDestroy(pme_event_finish));
   }

   if (op & rc_alloc) {
      g::q0 = acc_get_default_async();
      g::q1 = acc_async_sync;
      g::s0 = (cudaStream_t)acc_get_cuda_stream(g::q0);
      g::s1 = (cudaStream_t)acc_get_cuda_stream(g::q1);
      check_rt(cublasCreate(&g::h0)); // calls cudaMemcpy [sync] here
      check_rt(cublasCreate(&g::h1)); // calls cudaMemcpy [sync] here
      check_rt(cublasSetStream(g::h0, g::s0));
      check_rt(cublasSetStream(g::h1, g::s1));
      // set pointer mode for cublas dot kernels
      cublasPointerMode_t ptrflag = CUBLAS_POINTER_MODE_DEVICE;
      check_rt(cublasSetPointerMode(g::h0, ptrflag));
      check_rt(cublasSetPointerMode(g::h1, ptrflag));

      int nblock = get_grid_size(BLOCK_DIM);
      check_rt(cudaMallocHost(&pinned_buf, nblock * sizeof(double)));
      check_rt(cudaMalloc(&dptr_buf, nblock * sizeof(double)));

      use_pme_stream = false;
      if (pltfm_config & CU_PLTFM) {
         g::qpme = g::q1;
         g::spme = g::s1;
      } else {
         g::qpme = g::q0 + 1;
         g::spme = (cudaStream_t)acc_get_cuda_stream(g::qpme);
      }
      check_rt(cudaEventCreateWithFlags(&pme_event_start, cudaEventDisableTiming));
      check_rt(cudaEventCreateWithFlags(&pme_event_finish, cudaEventDisableTiming));

      check_rt(cudaProfilerStart());
   }
#else
   if (op & rc_alloc) {
      use_pme_stream = false;
   }
#endif
}
}
