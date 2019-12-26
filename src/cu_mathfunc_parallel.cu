#include "cudalib.h"
#include "error.h"
#include "mathfunc_parallel_cu.h"
#include <cublas_v2.h>
#include <cuda_runtime.h>


TINKER_NAMESPACE_BEGIN
namespace platform {
namespace cu {
template <>
void dotprod<float>(float* ans, const float* a, const float* b, int nelem,
                    int sync)
{
   cublasHandle_t hd = (sync ? h_cublas : h_cublas_nonblk);
   float alpha = 1, beta = 0;
   check_rt(cublasSgemm(hd, CUBLAS_OP_N, CUBLAS_OP_T, 1, 1, nelem, //
                        &alpha, a, 1, b, 1,                        //
                        &beta, ans, 1));
   if (sync)
      check_rt(cudaStreamSynchronize(nullptr));
}


template <>
void dotprod<double>(double* ans, const double* a, const double* b, int nelem,
                     int sync)
{
   cublasHandle_t hd = (sync ? h_cublas : h_cublas_nonblk);
   double alpha = 1, beta = 0;
   check_rt(cublasDgemm(hd, CUBLAS_OP_N, CUBLAS_OP_T, 1, 1, nelem, //
                        &alpha, a, 1, b, 1,                        //
                        &beta, ans, 1));
   if (sync)
      check_rt(cudaStreamSynchronize(nullptr));
}
}
}
TINKER_NAMESPACE_END
