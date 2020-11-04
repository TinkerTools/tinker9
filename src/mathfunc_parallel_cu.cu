#include "mathfunc_parallel_cu.h"
#include "syntax/cu/reduce.h"
#include "tool/cudalib.h"
#include "tool/deduce_ptr.h"
#include "tool/error.h"
#include "tool/gpu_card.h"
#include <cassert>
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <numeric>


namespace tinker {
namespace {
template <class T, class Op>
void reduce_to_dptr(const T* a, size_t nelem, cudaStream_t st)
{
   T* dptr = (T*)dptr_buf;
   int grid_siz1 = get_grid_size(BLOCK_DIM);
   int grid_siz2 = (nelem + BLOCK_DIM - 1) / BLOCK_DIM;
   int grid_size = std::min(grid_siz1, grid_siz2);
   reduce<T, BLOCK_DIM, Op><<<grid_size, BLOCK_DIM, 0, st>>>(dptr, a, nelem);
   reduce<T, BLOCK_DIM, Op><<<1, BLOCK_DIM, 0, st>>>(dptr, dptr, grid_size);
}


template <class T, class Op>
T reduce_general(const T* a, size_t nelem, int queue)
{
   cudaStream_t st = queue == g::q1 ? g::s1 : g::s0;
   T* dptr = (T*)dptr_buf;
   T* hptr = (T*)pinned_buf;
   reduce_to_dptr<T, Op>(a, nelem, st);
   check_rt(cudaMemcpyAsync(hptr, dptr, sizeof(T), cudaMemcpyDeviceToHost, st));
   // always wait
   check_rt(cudaStreamSynchronize(st));
   return *hptr;
}
}


template <class T>
T reduce_sum_cu(const T* a, size_t nelem, int queue)
{
   return reduce_general<T, OpPlus<T>>(a, nelem, queue);
}
template int reduce_sum_cu(const int*, size_t, int);
template float reduce_sum_cu(const float*, size_t, int);
template double reduce_sum_cu(const double*, size_t, int);
template unsigned long long reduce_sum_cu(const unsigned long long*, size_t,
                                          int);


template <class HT, size_t HN, class DPTR>
void reduce_sum2_cu(HT (&restrict h_ans)[HN], DPTR restrict a, size_t nelem,
                    int queue)
{
   typedef typename deduce_ptr<DPTR>::type CONST_DT;
   typedef typename std::remove_const<CONST_DT>::type T;
   static_assert(std::is_same<HT, T>::value, "");
   constexpr size_t N = deduce_ptr<DPTR>::n;
   static_assert(HN <= N, "");

   cudaStream_t st = queue == g::q1 ? g::s1 : g::s0;
   T(*dptr)[HN] = (T(*)[HN])dptr_buf;
   T* hptr = (T*)pinned_buf;
   int grid_siz1 = get_grid_size(BLOCK_DIM);
   grid_siz1 = grid_siz1 / HN; // limited by the output buffer
   int grid_siz2 = (nelem + BLOCK_DIM - 1) / BLOCK_DIM;
   int grid_size = std::min(grid_siz1, grid_siz2);
   reduce2<T, BLOCK_DIM, HN, N, OpPlus<T>>
      <<<grid_size, BLOCK_DIM, 0, st>>>(dptr, a, nelem);
   reduce2<T, BLOCK_DIM, HN, HN, OpPlus<T>>
      <<<1, BLOCK_DIM, 0, st>>>(dptr, dptr, grid_size);
   check_rt(cudaMemcpyAsync(hptr, (T*)dptr, HN * sizeof(HN),
                            cudaMemcpyDeviceToHost, st));
   // always wait
   check_rt(cudaStreamSynchronize(st));
   #pragma unroll
   for (int j = 0; j < HN; ++j)
      h_ans[j] = hptr[j];
}
template void reduce_sum2_cu(float (&)[6], float (*)[8], size_t, int);
template void reduce_sum2_cu(double (&)[6], double (*)[8], size_t, int);
template void reduce_sum2_cu(unsigned long long (&)[6],
                             unsigned long long (*)[8], size_t, int);


template <>
void dotprod_cu<float>(float* ans, const float* a, const float* b, size_t nelem,
                       int queue)
{
   bool dq = queue == g::q1;
   cublasHandle_t hd = (dq ? g::h1 : g::h0);
   check_rt(cublasSdot(hd, nelem, a, 1, b, 1, ans));
}


template <>
void dotprod_cu<double>(double* ans, const double* a, const double* b,
                        size_t nelem, int queue)
{
   bool dq = queue == g::q1;
   cublasHandle_t hd = (dq ? g::h1 : g::h0);
   check_rt(cublasDdot(hd, nelem, a, 1, b, 1, ans));
}


// cublas gemm does not run as fast here prior to cuda 10.1.
// Old code:
//
// #if CUDART_VERSION >= 10100 // >= 10.1
//    float alpha = 1, beta = 0;
//    check_rt(cublasSgemm(hd, CUBLAS_OP_N, CUBLAS_OP_T, 1, 1, nelem, //
//                         &alpha, a, 1, b, 1,                        //
//                         &beta, ans, 1));
// #else
//    check_rt(cublasSdot(hd, nelem, a, 1, b, 1, ans));
// #endif
}
