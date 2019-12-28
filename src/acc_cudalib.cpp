#include "acclib.h"
#include "cudalib.h"
#if TINKER_CUDART
#   include "error.h"
#   include "gpu_card.h"
#   include <openacc.h>
#endif


TINKER_NAMESPACE_BEGIN
#if TINKER_CUDART
cudaStream_t nonblk;
cublasHandle_t h_cublas;
cublasHandle_t h_cublas_nonblk;
real* pinned_real64;
real* dptr_real64;
#endif


int async_queue;
void cudalib_data(rc_op op)
{
#if TINKER_CUDART
   if (op & rc_dealloc) {
      async_queue = -42;
      nonblk = nullptr;
      check_rt(cublasDestroy(h_cublas));
      check_rt(cublasDestroy(h_cublas_nonblk));
      check_rt(cudaFreeHost(pinned_real64));
      check_rt(cudaFree(dptr_real64));
   }


   if (op & rc_alloc) {
      async_queue = acc_get_default_async();
      nonblk = (cudaStream_t)acc_get_cuda_stream(async_queue);
      check_rt(cublasCreate(&h_cublas));
      check_rt(cublasCreate(&h_cublas_nonblk));
      check_rt(cublasSetStream(h_cublas_nonblk, nonblk));
      int nblock = get_grid_size(BLOCK_DIM);
      check_rt(cudaMallocHost(&pinned_real64, nblock * sizeof(double)));
      check_rt(cudaMalloc(&dptr_real64, nblock * sizeof(double)));
   }
#endif
}


void wait_queue()
{
#if TINKER_CUDART
#pragma acc wait(async_queue)
#endif
}
TINKER_NAMESPACE_END
