#include "cudalib.h"
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


template <unsigned int B, class T, class Op>
__device__
void warp_reduce(volatile T* sdata, unsigned int t, Op op)
{
   // clang-format off
   if (B >= 64) sdata[t] = op(sdata[t], sdata[t + 32]);
   if (B >= 32) sdata[t] = op(sdata[t], sdata[t + 16]);
   if (B >= 16) sdata[t] = op(sdata[t], sdata[t + 8] );
   if (B >= 8)  sdata[t] = op(sdata[t], sdata[t + 4] );
   if (B >= 4)  sdata[t] = op(sdata[t], sdata[t + 2] );
   if (B >= 2)  sdata[t] = op(sdata[t], sdata[t + 1] );
   // clang-format on
}


template <unsigned int B, class T, class Op>
__global__
void reduce11(T* g_odata, const T* g_idata, size_t n, Op op = Op())
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
   if (B >= 128) { if (t < 64)  { sd[t] = op(sd[t], sd[t + 64] ); } SYNC(); }
   if (t < 32) warp_reduce<B, T, Op>(sd, t, op);
   if (t == 0) g_odata[blockIdx.x] = sd[0];
   // clang-format on
}


template <class T, class Op>
void reduce_to_dptr(const T* a, size_t nelem, int sync)
{
   cudaStream_t st = (sync ? nullptr : nonblk);
   T* dptr = (T*)dptr_real64;
   int grid_size = get_grid_size(BLOCK_DIM);
   reduce11<BLOCK_DIM, T, Op><<<grid_size, BLOCK_DIM, 0, st>>>(dptr, a, nelem);
   reduce11<BLOCK_DIM, T, Op><<<1, BLOCK_DIM, 0, st>>>(dptr, dptr, grid_size);
}


template <class T, class Op>
T reduce_general(const T* a, size_t nelem, int sync)
{
   cudaStream_t st = (sync ? nullptr : nonblk);
   T* dptr = (T*)dptr_real64;
   T* hptr = (T*)pinned_real64;
   int grid_size = get_grid_size(BLOCK_DIM);
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
