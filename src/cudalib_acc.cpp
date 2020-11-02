#include "glob.accasync.h"
#include "tool/cudalib.h"
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


      sync_queue = -42;
      async_queue = -42;
      nonblk = nullptr;
      check_rt(cublasDestroy(h_cublas));
      check_rt(cublasDestroy(h_cublas_nonblk));
      check_rt(cudaFreeHost(pinned_buf));
      check_rt(cudaFree(dptr_buf));


      use_pme_stream = false;
      pme_stream = nullptr;
      pme_queue = -42;
      check_rt(cudaEventDestroy(pme_event_finish));
   }


   if (op & rc_alloc) {
      sync_queue = acc_async_sync;
      async_queue = acc_get_default_async();
      nonblk = (cudaStream_t)acc_get_cuda_stream(async_queue);
      check_rt(cublasCreate(&h_cublas));        // calls cudaMemcpy [sync] here
      check_rt(cublasCreate(&h_cublas_nonblk)); // calls cudaMemcpy [sync] here
      check_rt(cublasSetStream(h_cublas_nonblk, nonblk));
      // set pointer mode for cublas dot kernels
      cublasPointerMode_t ptrflag = CUBLAS_POINTER_MODE_DEVICE;
      check_rt(cublasSetPointerMode(h_cublas, ptrflag));
      check_rt(cublasSetPointerMode(h_cublas_nonblk, ptrflag));


      int nblock = get_grid_size(BLOCK_DIM);
      check_rt(cudaMallocHost(&pinned_buf, nblock * sizeof(double)));
      check_rt(cudaMalloc(&dptr_buf, nblock * sizeof(double)));


      use_pme_stream = false;
      pme_queue = async_queue + 1;
      pme_stream = (cudaStream_t)acc_get_cuda_stream(pme_queue);
      check_rt(
         cudaEventCreateWithFlags(&pme_event_finish, cudaEventDisableTiming));


      check_rt(cudaProfilerStart());
   }
#else
   if (op & rc_alloc) {
      use_pme_stream = false;
      pme_queue = -42;
   }
#endif
}
}
