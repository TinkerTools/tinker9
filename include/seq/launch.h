#pragma once
#include "ff/energy.h"
#include "tool/cudalib.h"
#include "tool/gpucard.h"

namespace tinker {
/// \addtogroup cuda_syntax
///\{

/// Launch a non-blocking CUDA kernel.
///
/// \param st  CUDA stream.
/// \param sh  Dynamic shared memory (bytes per block).
/// \param bs  Thread block size.
/// \param np  Hint to determine the grid dimension.
///            This is often the loop limit.
/// \param k   CUDA \c __global__ kernel.
/// \param a   Arguments of the kernel.
///
/// \see launch_k3b
template <class K, class... Ts>
void launch_k3s(cudaStream_t st, size_t sh, int bs, int np, K k, Ts&&... a)
{
   int gs = (np + bs - 1) / bs;
   k<<<gs, bs, sh, st>>>(std::forward<Ts>(a)...);
}

/// Launch a non-blocking CUDA kernel.
/// Energy, virial, and/or count buffers are involved.
///
/// \param st  CUDA stream.
/// \param sh  Dynamic shared memory (bytes per block).
/// \param bs  Thread block size.
/// \param np  Hint to determine the grid dimension.
///            This is often the loop limit.
/// \param k   CUDA \c __global__ kernel.
/// \param a   Arguments of the kernel.
///
/// \see launch_k3s
template <class K, class... Ts>
void launch_k3b(cudaStream_t st, size_t sh, int bs, int np, K k, Ts&&... a)
{
   np = std::min(np, nelem_buffer);
   int gs = (np + bs - 1) / bs;
   k<<<gs, bs, sh, st>>>(std::forward<Ts>(a)...);
}

/// Launch a non-blocking CUDA kernel with 0 dynamic shared memory.
///
/// \see launch_k3s
template <class K, class... Ts>
void launch_k2s(cudaStream_t st, int bs, int np, K k, Ts&&... a)
{
   const size_t sh = 0;
   launch_k3s(st, sh, bs, np, k, std::forward<Ts>(a)...);
}

/// Launch a non-blocking CUDA kernel with 0 dynamic shared memory.
/// Energy, virial, and/or count buffers are involved.
///
/// \see launch_k3b
template <class K, class... Ts>
void launch_k2b(cudaStream_t st, int bs, int np, K k, Ts&&... a)
{
   const size_t sh = 0;
   launch_k3b(st, sh, bs, np, k, std::forward<Ts>(a)...);
}

/// Launch a non-blocking CUDA kernel with 0 dynamic shared memory
/// and the default block size (#BLOCK_DIM).
///
/// \see BLOCK_DIM, launch_k2s, launch_k3s
template <class K, class... Ts>
void launch_k1s(cudaStream_t st, int np, K k, Ts&&... a)
{
   const int bs = BLOCK_DIM;
   launch_k2s(st, bs, np, k, std::forward<Ts>(a)...);
}

/// Launch a non-blocking CUDA kernel with 0 dynamic shared memory
/// and the default block size (#BLOCK_DIM).
/// Energy, virial, and/or count buffers are involved.
///
/// \see BLOCK_DIM, launch_k2b, launch_k3b
template <class K, class... Ts>
void launch_k1b(cudaStream_t st, int np, K k, Ts&&... a)
{
   const int bs = BLOCK_DIM;
   launch_k2b(st, bs, np, k, std::forward<Ts>(a)...);
}

/// \def ITHREAD
/// \def STRIDE
#define ITHREAD threadIdx.x + blockIdx.x* blockDim.x
#define STRIDE  blockDim.x* gridDim.x

/// \}
}
