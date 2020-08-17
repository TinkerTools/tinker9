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


      async_queue = -42;
      nonblk = nullptr;
      check_rt(cudaStreamDestroy(stream2));
      check_rt(cublasDestroy(h_cublas));
      check_rt(cublasDestroy(h_cublas_nonblk));
      check_rt(cudaFreeHost(pinned_buf));
      check_rt(cudaFree(dptr_buf));


      cudaEventDestroy(stream2_begin_event);
      cudaEventDestroy(stream2_end_event);
      use_stream2 = false;
   }


   if (op & rc_alloc) {
      async_queue = acc_get_default_async();
      nonblk = (cudaStream_t)acc_get_cuda_stream(async_queue);
      check_rt(cudaStreamCreateWithFlags(&stream2, cudaStreamNonBlocking));
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


      cudaEventCreateWithFlags(&stream2_begin_event, cudaEventDisableTiming);
      cudaEventCreateWithFlags(&stream2_end_event, cudaEventDisableTiming);
      use_stream2 = false;


      check_rt(cudaProfilerStart());
   }
#else
   (void)op;
#endif
}


void stream2_sync()
{
#if TINKER_CUDART
   if (use_stream2) {
      check_rt(cudaEventRecord(stream2_end_event, stream2));
      check_rt(cudaStreamWaitEvent(nonblk, stream2_end_event, 0));
   }
#endif
}


void stream2_begin()
{
#if TINKER_CUDART
   if (use_stream2) {
      // Record `stream2_begin_event` when other kernels on `nonblk` have ended.
      check_rt(cudaEventRecord(stream2_begin_event, nonblk));
      // `stream2` will wait until `stream2_begin_event` is recorded.
      check_rt(cudaStreamWaitEvent(stream2, stream2_begin_event, 0));
   }
#endif
}
}
