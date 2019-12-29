#include "cudalib.h"
#include "deduce_ptr.h"
#include "error.h"
#include "gpu_card.h"
#include "mathfunc_parallel_cu.h"
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <numeric>


TINKER_NAMESPACE_BEGIN
namespace platform {
namespace cu {
template <class T>
struct OpPlus
{
   __device__
   T init() const
   {
      return 0;
   }


   __device__
   T operator()(T a, T b) const
   {
      return a + b;
   }
};


template <class T>
struct OpLogicOr
{
   __device__
   T init() const
   {
      return false;
   }


   __device__
   T operator()(T a, T b) const
   {
      return a || b;
   }
};


template <class T, unsigned int B, class Op>
__device__
void warp_reduce(volatile T* sd, unsigned int t, Op op)
{
   // clang-format off
   if (B >= 64) sd[t] = op(sd[t], sd[t + 32]);
   if (B >= 32) sd[t] = op(sd[t], sd[t + 16]);
   if (B >= 16) sd[t] = op(sd[t], sd[t + 8 ]);
   if (B >= 8)  sd[t] = op(sd[t], sd[t + 4 ]);
   if (B >= 4)  sd[t] = op(sd[t], sd[t + 2 ]);
   if (B >= 2)  sd[t] = op(sd[t], sd[t + 1 ]);
   // clang-format on
}


template <class T, unsigned int HN, unsigned int B, class Op>
__device__
void warp_reduce2(volatile T (*sd)[B], unsigned int t, Op op)
{
#define RDC2_ADD(x)                                                            \
   _Pragma("unroll") for (int j = 0; j < HN; ++j) sd[j][t] =                   \
      op(sd[j][t], sd[j][t + x])
   // clang-format off
   if (B >= 64) RDC2_ADD(32);
   if (B >= 32) RDC2_ADD(16);
   if (B >= 16) RDC2_ADD(8 );
   if (B >= 8)  RDC2_ADD(4 );
   if (B >= 4)  RDC2_ADD(2 );
   if (B >= 2)  RDC2_ADD(1 );
   // clang-format on
}


template <class T, unsigned int B, class Op>
__global__
void reduce(T* g_odata, const T* g_idata, size_t n, Op op = Op())
{
   __shared__ T sd[B];
   unsigned int t = threadIdx.x;
   sd[t] = op.init();
   for (int i = t + blockIdx.x * B; i < n; i += B * gridDim.x) {
      sd[t] = op(sd[t], g_idata[i]);
   }
   __syncthreads();


#define SYNC __syncthreads
   // clang-format off
   if (B >= 512) { if (t < 256) { sd[t] = op(sd[t], sd[t + 256]); } SYNC(); }
   if (B >= 256) { if (t < 128) { sd[t] = op(sd[t], sd[t + 128]); } SYNC(); }
   if (B >= 128) { if (t < 64 ) { sd[t] = op(sd[t], sd[t + 64 ]); } SYNC(); }
   if (t < 32  ) warp_reduce<T, B, Op>(sd, t, op);
   // clang-format on
   if (t == 0)
      g_odata[blockIdx.x] = sd[0];
}


template <class T, unsigned int B, unsigned int HN, size_t N, class Op>
__global__
void reduce2(T (*g_odata)[HN], const T (*g_idata)[N], size_t n, Op op = Op())
{
   __shared__ T sd[HN][B];
   unsigned int t = threadIdx.x;
   #pragma unroll
   for (int j = 0; j < HN; ++j)
      sd[j][t] = 0;
   for (int i = t + blockIdx.x * B; i < n; i += B * gridDim.x) {
      #pragma unroll
      for (int j = 0; j < HN; ++j)
         sd[j][t] = op(sd[j][t], g_idata[i][j]);
   }
   __syncthreads();


   // clang-format off
   if (B >= 512) { if (t < 256) { RDC2_ADD(256); } SYNC(); }
   if (B >= 256) { if (t < 128) { RDC2_ADD(128); } SYNC(); }
   if (B >= 128) { if (t < 64 ) { RDC2_ADD(64 ); } SYNC(); }
   if (t < 32  ) warp_reduce2<T, HN, B, Op>(sd, t, op);
   // clang-format on
   if (t == 0)
      #pragma unroll
      for (int j = 0; j < HN; ++j)
         g_odata[blockIdx.x][j] = sd[j][0];
}


template <class T, class Op>
void reduce_to_dptr(const T* a, size_t nelem, int sync)
{
   cudaStream_t st = (sync ? nullptr : nonblk);
   T* dptr = (T*)dptr_real64;
   int grid_siz1 = get_grid_size(BLOCK_DIM);
   int grid_siz2 = (nelem + BLOCK_DIM - 1) / BLOCK_DIM;
   int grid_size = std::min(grid_siz1, grid_siz2);
   reduce<T, BLOCK_DIM, Op><<<grid_size, BLOCK_DIM, 0, st>>>(dptr, a, nelem);
   reduce<T, BLOCK_DIM, Op><<<1, BLOCK_DIM, 0, st>>>(dptr, dptr, grid_size);
}


template <class T, class Op>
T reduce_general(const T* a, size_t nelem, int sync)
{
   cudaStream_t st = (sync ? nullptr : nonblk);
   T* dptr = (T*)dptr_real64;
   T* hptr = (T*)pinned_real64;
   reduce_to_dptr<T, Op>(a, nelem, sync);
   check_rt(cudaMemcpyAsync(hptr, dptr, sizeof(T), cudaMemcpyDeviceToHost, st));
   check_rt(cudaStreamSynchronize(st));
   return *hptr;
}


template <class T>
T reduce_sum(const T* a, size_t nelem, int sync)
{
   return reduce_general<T, OpPlus<T>>(a, nelem, sync);
}
template int reduce_sum(const int*, size_t, int);
template float reduce_sum(const float*, size_t, int);
template double reduce_sum(const double*, size_t, int);
template unsigned long long reduce_sum(const unsigned long long*, size_t, int);


template <class HT, size_t HN, class DPTR>
void reduce_sum2(HT (&restrict h_ans)[HN], DPTR restrict a, size_t nelem,
                 int sync)
{
   typedef typename deduce_ptr<DPTR>::type CONST_DT;
   typedef typename std::remove_const<CONST_DT>::type T;
   static_assert(std::is_same<HT, T>::value, "");
   constexpr size_t N = deduce_ptr<DPTR>::n;
   static_assert(HN <= N, "");


   cudaStream_t st = (sync ? nullptr : nonblk);
   T(*dptr)[HN] = (T(*)[HN])dptr_real64;
   T* hptr = (T*)pinned_real64;
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
   check_rt(cudaStreamSynchronize(st));
   #pragma unroll
   for (int j = 0; j < HN; ++j)
      h_ans[j] = hptr[j];
}
template void reduce_sum2(float (&)[6], float (*)[8], size_t, int);
template void reduce_sum2(double (&)[6], double (*)[8], size_t, int);
template void reduce_sum2(unsigned long long (&)[6], unsigned long long (*)[8],
                          size_t, int);


template <class T>
T reduce_logic_or(const T* a, size_t nelem, int sync)
{
   return reduce_general<T, OpLogicOr<T>>(a, nelem, sync);
}
template int reduce_logic_or(const int*, size_t, int);


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
