#include "macro.h"


#if TINKER_CUDART && defined(_OPENACC)
#   include "cudalib.h"
#   include "dev_array.h"
#   include "io_print.h"
#   include <openacc.h>


TINKER_NAMESPACE_BEGIN
cudaStream_t async_acc;
cudaStream_t nonblk;


cublasHandle_t h_cublas;
cublasHandle_t h_cublas_nonblk;
real* pinned_real64;
real* dptr_real64;


void cudalib_data(rc_op op)
{
   if (op & rc_dealloc) {
      async_acc = nullptr;
      nonblk = nullptr;


      check_rt(cublasDestroy(h_cublas));
      check_rt(cublasDestroy(h_cublas_nonblk));
      check_rt(cudaFreeHost(pinned_real64));
      device_array::deallocate(dptr_real64);
   }


   if (op & rc_alloc) {
      int handle;
      // const char* fmt = " {:30s}{:>12}\n";


      // handle = acc_async_sync;
      // print(stdout, fmt, "acc async sync:", handle);
      handle = acc_get_default_async();
      // print(stdout, fmt, "acc default async:", handle);


      async_acc = (cudaStream_t)acc_get_cuda_stream(acc_async_sync);
      // print(stdout, fmt, "cuda stream1 (async_acc):", (void*)async_acc);
      nonblk = (cudaStream_t)acc_get_cuda_stream(handle);
      // print(stdout, fmt, "cuda stream2 (nonblk):", (void*)nonblk);


      check_rt(cublasCreate(&h_cublas));
      check_rt(cublasSetStream(h_cublas, async_acc));
      check_rt(cublasCreate(&h_cublas_nonblk));
      check_rt(cublasSetStream(h_cublas_nonblk, nonblk));
      check_rt(cudaMallocHost(&pinned_real64, 64 * sizeof(real)));
      device_array::allocate(64, &dptr_real64);
   }
}
TINKER_NAMESPACE_END
#endif
