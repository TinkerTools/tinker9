#pragma once
#include "macro.h"


TINKER_NAMESPACE_BEGIN
namespace pltfm_cu {
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
   // clang-format off
   if (B >= 64) _Pragma("unroll") for (int j = 0; j < HN; ++j) sd[j][t] = op(sd[j][t], sd[j][t + 32]);
   if (B >= 32) _Pragma("unroll") for (int j = 0; j < HN; ++j) sd[j][t] = op(sd[j][t], sd[j][t + 16]);
   if (B >= 16) _Pragma("unroll") for (int j = 0; j < HN; ++j) sd[j][t] = op(sd[j][t], sd[j][t + 8 ]);
   if (B >= 8)  _Pragma("unroll") for (int j = 0; j < HN; ++j) sd[j][t] = op(sd[j][t], sd[j][t + 4 ]);
   if (B >= 4)  _Pragma("unroll") for (int j = 0; j < HN; ++j) sd[j][t] = op(sd[j][t], sd[j][t + 2 ]);
   if (B >= 2)  _Pragma("unroll") for (int j = 0; j < HN; ++j) sd[j][t] = op(sd[j][t], sd[j][t + 1 ]);
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


   // clang-format off
   if (B >= 512) { if (t < 256) { sd[t] = op(sd[t], sd[t + 256]); } __syncthreads(); }
   if (B >= 256) { if (t < 128) { sd[t] = op(sd[t], sd[t + 128]); } __syncthreads(); }
   if (B >= 128) { if (t < 64 ) { sd[t] = op(sd[t], sd[t + 64 ]); } __syncthreads(); }
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
   if (B >= 512) { if (t < 256) { _Pragma("unroll") for (int j = 0; j < HN; ++j) sd[j][t] = op(sd[j][t], sd[j][t + 256]); } __syncthreads(); }
   if (B >= 256) { if (t < 128) { _Pragma("unroll") for (int j = 0; j < HN; ++j) sd[j][t] = op(sd[j][t], sd[j][t + 128]); } __syncthreads(); }
   if (B >= 128) { if (t < 64 ) { _Pragma("unroll") for (int j = 0; j < HN; ++j) sd[j][t] = op(sd[j][t], sd[j][t + 64 ]); } __syncthreads(); }
   if (t < 32  ) warp_reduce2<T, HN, B, Op>(sd, t, op);
   // clang-format on
   if (t == 0)
      #pragma unroll
      for (int j = 0; j < HN; ++j)
         g_odata[blockIdx.x][j] = sd[j][0];
}
}
TINKER_NAMESPACE_END
