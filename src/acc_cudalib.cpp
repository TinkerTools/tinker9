#include "cudalib.h"
#if TINKER_CUDART
#   include "error.h"
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


void cudalib_data(rc_op op)
{
#if TINKER_CUDART
   if (op & rc_dealloc) {
      nonblk = nullptr;
      check_rt(cublasDestroy(h_cublas));
      check_rt(cublasDestroy(h_cublas_nonblk));
      check_rt(cudaFreeHost(pinned_real64));
      check_rt(cudaFree(dptr_real64));
   }


   if (op & rc_alloc) {
      int handle = acc_get_default_async();
      nonblk = (cudaStream_t)acc_get_cuda_stream(handle);

      check_rt(cublasCreate(&h_cublas));
      check_rt(cublasCreate(&h_cublas_nonblk));
      check_rt(cublasSetStream(h_cublas_nonblk, nonblk));
      check_rt(cudaMallocHost(&pinned_real64, 64 * sizeof(real)));
      check_rt(cudaMalloc(&dptr_real64, 64 * sizeof(real)));
   }
#endif
}
TINKER_NAMESPACE_END
