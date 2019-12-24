#include "cudalib.h"
#include "cudalib2.h"


TINKER_NAMESPACE_BEGIN
template <class A, class T1, class T2>
void dotprod(cudaStream_t st, A* ans, const T1* a, const T2* b, int nelem)
{
   cublasHandle_t hd = (st == nullptr ? h_cublas : h_cublas_nonblk);
   if CONSTEXPR (std::is_floating_point<A>::value) {
      if CONSTEXPR (sizeof(real) == sizeof(float)) {
         float alpha = 1;
         float beta = 0;
         check_rt(cublasSgemm(hd, CUBLAS_OP_N, CUBLAS_OP_T, 1, 1, nelem, //
                              &alpha, (const float*)a, 1,                //
                              (const float*)b, 1,                        //
                              &beta, (float*)ans, 1));
         if (st == nullptr)
            check_rt(cudaStreamSynchronize(nullptr));
         return;
      } else if CONSTEXPR (sizeof(real) == sizeof(double)) {
         double alpha = 1;
         double beta = 0;
         check_rt(cublasDgemm(hd, CUBLAS_OP_N, CUBLAS_OP_T, 1, 1, nelem, //
                              &alpha, (const double*)a, 1,               //
                              (const double*)b, 1,                       //
                              &beta, (double*)ans, 1));
         if (st == nullptr)
            check_rt(cudaStreamSynchronize(nullptr));
         return;
      }
   }


   assert(false);
}


template void dotprod(cudaStream_t, float*, const float (*)[3],
                      const float (*)[3], int);
template void dotprod(cudaStream_t, double*, const double (*)[3],
                      const double (*)[3], int);
TINKER_NAMESPACE_END
