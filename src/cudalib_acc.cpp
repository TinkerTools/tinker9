#include "acclib.h"
#include "cudalib.h"
#include "wait_queue.h"
#if TINKER_CUDART
#   include "error.h"
#   include "gpu_card.h"
#   include <cuda_profiler_api.h>
#   include <openacc.h>
#endif


namespace tinker {
#if TINKER_CUDART
cudaStream_t nonblk;
cublasHandle_t h_cublas;
cublasHandle_t h_cublas_nonblk;
void* pinned_buf;
void* dptr_buf;
#endif


int async_queue;
void cudalib_data(rc_op op)
{
#if TINKER_CUDART
   if (op & rc_dealloc) {
      check_rt(cudaProfilerStop());


      async_queue = -42;
      nonblk = nullptr;
      check_rt(cublasDestroy(h_cublas));
      check_rt(cublasDestroy(h_cublas_nonblk));
      check_rt(cudaFreeHost(pinned_buf));
      check_rt(cudaFree(dptr_buf));
   }


   if (op & rc_alloc) {
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


      check_rt(cudaProfilerStart());
   }
#else
   (void)op;
#endif
}


void wait_queue(DMFlag flag)
{
#if TINKER_CUDART
   if ((flag & DMFlag::WAIT) && !(flag & DMFlag::DEFAULT_Q)) {
      #pragma acc wait(async_queue)
   }
#else
   (void)flag;
#endif
}
}
